import os
import subprocess


def run(dir):
    for filename in os.listdir(dir):
        subprocess.run(["bash", os.path.join(dir, filename)])