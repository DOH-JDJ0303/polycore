from pathlib import Path
import subprocess, sys

def test_cli_help():
    r = subprocess.run([sys.executable, "-m", "polycore", "--version"], capture_output=True)
    assert r.returncode == 0
