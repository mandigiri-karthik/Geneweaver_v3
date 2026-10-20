"""Annotate the frozen Phase 0 reference genome with the Bakta web API.

Usage:
    python3 run_bakta.py <fasta_file> <output_folder>

Steps: create job -> upload files -> start job -> print logs until done -> download results.
Needs only Python 3 and curl (both already on macOS).

Only the result files used for scoring are downloaded (no plots, EMBL or GenBank).
The job ID is written to job_id.txt; the job secret is never written to disk.
"""
import json
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://api.bakta.computational.bio/api/v1"
CHECK_EVERY_SECONDS = 30

# Bakta result types to download (keys of the job/result reply)
WANTED_RESULTS = {"GFF3", "TSV", "TSVHypothetical", "TSVInference", "FAA", "JSON", "TXTLogs"}


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


def get_version():
    """Return the server's Bakta tool and database versions, if the API reports them."""
    try:
        with urllib.request.urlopen(f"{API}/version") as response:
            return json.loads(response.read().decode())
    except Exception as error:  # the endpoint is informational only
        return {"error": str(error)}


def create_job(name):
    """Create a new job. Returns the job (jobID + secret) and the upload links."""
    reply = post_json("job/init", {"name": name, "repliconTableType": "TSV"})
    print("Job ID:", reply["job"]["jobID"])
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
    print("Job started with config:", json.dumps(config))


def get_status(job):
    """Return the job status, e.g. INIT, RUNNING, SUCCESSFUL or ERROR.

    Right after job/start the server can briefly return an empty job list,
    so that case is reported as PENDING instead of failing.
    """
    reply = post_json("job/list", {"jobs": [job]})
    if not reply or not reply.get("jobs"):
        return "PENDING"
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
        print(f"[{time.strftime('%H:%M:%S')}] status: {status}", flush=True)
        time.sleep(CHECK_EVERY_SECONDS)


def download_results(job, output_folder):
    """Download the wanted result files into the output folder."""
    reply = post_json("job/result", job)
    for file_type, link in reply["ResultFiles"].items():
        if file_type not in WANTED_RESULTS:
            continue
        file_name = link.split("?")[0].split("/")[-1]
        urllib.request.urlretrieve(link, Path(output_folder) / file_name)
        print(f"Downloaded {file_type}: {file_name}")


def main():
    fasta_path = Path(sys.argv[1])
    output_folder = Path(sys.argv[2])
    output_folder.mkdir(parents=True, exist_ok=True)

    print("Started:", time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    print("Server version:", json.dumps(get_version()))

    reply = create_job(fasta_path.stem)
    job = reply["job"]
    (output_folder / "job_id.txt").write_text(job["jobID"] + "\n")

    replicon_table = output_folder / "replicons.tsv"
    make_replicon_table(fasta_path, replicon_table)
    upload_file(fasta_path, reply["uploadLinkFasta"])
    upload_file(replicon_table, reply["uploadLinkReplicons"])

    start_job(job)
    status = watch_job(job)
    print("Final status:", status)

    if status == "SUCCESSFUL":
        download_results(job, output_folder)
    print("Finished:", time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))


if __name__ == "__main__":
    main()
