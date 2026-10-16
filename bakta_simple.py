"""Annotate a genome with the Bakta web API and print the server logs live.

Usage:
    python3 bakta_simple.py <fasta_file> <output_folder>

Steps: create job -> upload files -> start job -> print logs until done -> download results.
Needs only Python 3 and curl (both already on macOS).
"""
import json
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://api.bakta.computational.bio/api/v1"
CHECK_EVERY_SECONDS = 5


def post_json(endpoint, body):
    """Send a JSON POST request to the Bakta API and return the JSON reply."""
    request = urllib.request.Request(
        f"{API}/{endpoint}",
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request) as response:
        text = response.read().decode()
    return json.loads(text) if text else None


def create_job(name):
    """Create a new job. Returns the job (jobID + secret) and the upload links."""
    reply = post_json("job/init", {"name": name, "repliconTableType": "TSV"})
    print("Job ID:", reply["job"]["jobID"])
    print("Secret:", reply["job"]["secret"])
    return reply


def upload_file(file_path, upload_link):
    """Upload a file to the link Bakta gave us.

    curl -T is used because it sends no Content-Type header,
    which the Bakta storage server requires.
    """
    subprocess.run(["curl", "-sS", "--fail", "-T", str(file_path), upload_link], check=True)
    print("Uploaded", file_path)


def make_replicon_table(fasta_path, table_path):
    """Write the replicon table Bakta expects: one line per sequence in the FASTA file.

    Every sequence is marked as a circular chromosome, which fits a complete bacterial genome.
    """
    lines = []
    for line in Path(fasta_path).read_text().splitlines():
        if line.startswith(">"):
            sequence_id = line[1:].split()[0]
            lines.append(f"{sequence_id}\t?\tchromosome\tcircular\t?\n")
    Path(table_path).write_text("".join(lines))


def start_job(job):
    """Tell Bakta to start annotating."""
    config = {
        "completeGenome": True,
        "translationTable": 11,  # genetic code for bacteria
        "keepContigHeaders": True,
        "minContigLength": 1,
        "compliant": False,
    }
    post_json("job/start", {"config": config, "job": job})
    print("Job started")


def get_status(job):
    """Return the job status, e.g. INIT, RUNNING, SUCCESSFUL or ERROR."""
    reply = post_json("job/list", {"jobs": [job]})
    return reply["jobs"][0]["jobStatus"]


def get_logs(job):
    """Return the full Bakta log text of the job."""
    query = urllib.parse.urlencode(job)
    with urllib.request.urlopen(f"{API}/job/logs?{query}") as response:
        return response.read().decode()


def watch_job(job):
    """Print new log lines as they appear, until the job finishes. Returns the final status."""
    lines_printed = 0
    while True:
        status = get_status(job)
        log_lines = get_logs(job).splitlines()

        for line in log_lines[lines_printed:]:
            print(line)
        lines_printed = len(log_lines)

        if status in ("SUCCESSFUL", "ERROR"):
            return status
        time.sleep(CHECK_EVERY_SECONDS)


def download_results(job, output_folder):
    """Download every result file into the output folder."""
    reply = post_json("job/result", job)
    for file_type, link in reply["ResultFiles"].items():
        file_name = link.split("?")[0].split("/")[-1]
        urllib.request.urlretrieve(link, Path(output_folder) / file_name)
        print(f"Downloaded {file_type}: {file_name}")


def main():
    fasta_path = Path(sys.argv[1])
    output_folder = Path(sys.argv[2])
    output_folder.mkdir(parents=True, exist_ok=True)

    reply = create_job(fasta_path.stem)
    job = reply["job"]

    replicon_table = output_folder / "replicons.tsv"
    make_replicon_table(fasta_path, replicon_table)
    upload_file(fasta_path, reply["uploadLinkFasta"])
    upload_file(replicon_table, reply["uploadLinkReplicons"])

    start_job(job)
    status = watch_job(job)
    print("Final status:", status)

    if status == "SUCCESSFUL":
        download_results(job, output_folder)


if __name__ == "__main__":
    main()
