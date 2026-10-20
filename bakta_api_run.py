"""Run Bakta genome annotation via the public Bakta web API.

Workflow (https://bakta.readthedocs.io/en/latest/api/full.html):
init -> upload FASTA -> start -> poll /job/list -> /job/result -> download files.
Uses only the Python standard library.
"""
import json
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

API = "https://api.bakta.computational.bio/api/v1"
FASTA = Path(sys.argv[1])
OUT_DIR = Path(sys.argv[2])
POLL_SECONDS = 30


def post(endpoint, body):
    req = urllib.request.Request(
        f"{API}/{endpoint}",
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as resp:
        raw = resp.read()
        return json.loads(raw) if raw else {}


def upload(path, url):
    # curl -T sends no Content-Type; the presigned URL rejects one
    code = subprocess.run(
        ["curl", "-sS", "-o", "/dev/null", "-w", "%{http_code}", "-T", str(path), url],
        check=True, capture_output=True, text=True,
    ).stdout
    if code != "200":
        sys.exit(f"Upload of {path.name} failed: HTTP {code}")
    print(f"Uploaded {path.name} (HTTP {code})")


def logs(job):
    query = urllib.parse.urlencode(job)
    with urllib.request.urlopen(f"{API}/job/logs?{query}", timeout=120) as resp:
        return resp.read().decode()


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    creds_file = OUT_DIR / "job.json"

    if creds_file.exists():
        job = json.loads(creds_file.read_text())
        print(f"Resuming existing job {job['jobID']}")
    else:
        # 1. Init
        init = post("job/init", {"name": FASTA.stem, "repliconTableType": "TSV"})
        job = init["job"]
        print(f"Initialized job {job['jobID']}")

        # 2. Upload FASTA and replicon table (Bakta errors on a missing table when repliconTableType is set)
        replicons = OUT_DIR / "replicons.tsv"
        replicons.write_text("".join(
            f"{line[1:].split()[0]}\t?\tchromosome\tcircular\t?\n"
            for line in FASTA.read_text().splitlines() if line.startswith(">")
        ))
        upload(FASTA, init["uploadLinkFasta"])
        upload(replicons, init["uploadLinkReplicons"])

        # 3. Start
        # The docs example only shows completeGenome, but the API requires all of these fields
        config = {
            "completeGenome": True,
            "translationTable": 11,  # bacterial genetic code
            "keepContigHeaders": True,
            "minContigLength": 1,
            "compliant": False,
        }
        post("job/start", {"config": config, "job": job})
        creds_file.write_text(json.dumps(job))
        print("Job started")

    # 4. Poll
    while True:
        status = post("job/list", {"jobs": [job]})
        if status.get("failedJobs"):
            sys.exit(f"Job lookup failed: {status['failedJobs']}")
        state = status["jobs"][0]["jobStatus"]
        print(f"[{time.strftime('%H:%M:%S')}] status: {state}", flush=True)
        if state == "SUCCESSFUL":
            break
        if state == "ERROR":
            print(logs(job))
            sys.exit("Bakta job ended with ERROR")
        time.sleep(POLL_SECONDS)

    # 5. Results
    result = post("job/result", job)
    (OUT_DIR / "result.json").write_text(json.dumps(result, indent=2))
    for kind, url in result["ResultFiles"].items():
        name = url.split("?")[0].rsplit("/", 1)[-1]
        urllib.request.urlretrieve(url, OUT_DIR / name)
        print(f"Downloaded {kind:16s} -> {name}")

    print("\n" + logs(job))

    # Summary from the annotation TSV
    tsv = next(OUT_DIR.glob("*.tsv"), None)
    for f in OUT_DIR.glob("*.tsv"):
        if not any(s in f.name for s in ("hypothetical", "inference")):
            tsv = f
    types = Counter()
    for line in tsv.read_text().splitlines():
        if line and not line.startswith("#"):
            types[line.split("\t")[1]] += 1
    print("\nFeature counts:")
    for t, n in types.most_common():
        print(f"  {t:12s} {n}")


if __name__ == "__main__":
    main()
