"""Entry point: runs the CAM agent. Same as `uv run scripts/run_cam.py ...`.

    uv run main.py --genome data/raw/lsakei_23k.fna --ref data/raw/lsakei_23k.gff3
"""
import runpy
import sys
from pathlib import Path

if __name__ == "__main__":
    sys.argv[0] = "run_cam.py"
    runpy.run_path(str(Path(__file__).parent / "scripts" / "run_cam.py"), run_name="__main__")
