#Run this script to run all scripts and save all results

import subprocess
import sys
from pathlib import Path

directory = Path(__file__).parent

scripts = ['advection.py', 
                'stability_test.py',
                'resolutionErrors.py']

for script in scripts:
    print(f"\nRunning {script}...")
    subprocess.run(
        [sys.executable, str(directory / script)],
        check=True
    )

print("\nAll scripts completed!")