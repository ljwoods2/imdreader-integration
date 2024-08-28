import imdreader
import MDAnalysis as mda
import pytest
import subprocess
import time
from pathlib import Path
import os
import logging

logger = logging.getLogger(__name__)
file_handler = logging.FileHandler("test.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)


### Utils, move to different file


def recvuntil(file_path, target_line, timeout):
    """
    Read from the file until the target line is found or the timeout occurs.

    Args:
        file_path (str): The path to the file to read from.
        target_line (str): The line to wait for.
        timeout (float): The timeout in seconds.

    Returns:
        str: The line containing the target line.

    Raises:
        TimeoutError: If the target line is not found within the timeout period.
    """
    end_time = time.time() + timeout
    with open(file_path, "r") as f:
        # Move to the end of the file initially
        f.seek(0, os.SEEK_END)
        while time.time() < end_time:
            time.sleep(0.1)  # Small delay to avoid busy-waiting
            # Read any new lines added to the file
            line = f.readline()
            if line:
                logger.debug(f"Read line: {line.strip()}")
                if target_line in line:
                    return line
            else:
                time.sleep(0.1)  # Delay if no new line is found
    raise TimeoutError(
        f"Timeout after {timeout} seconds waiting for '{target_line}'"
    )


def run_sim_and_wait(command, match_string, input_file=None):
    with open("sim_output.log", "w") as f:
        if input_file is None:
            p = subprocess.Popen(
                command,
                stdout=f,
                stderr=f,
                text=True,
                bufsize=1,
            )
        else:
            with open(input_file, "r") as input_f:
                p = subprocess.Popen(
                    command,
                    stdin=input_f,
                    stdout=f,
                    stderr=f,
                    text=True,
                    bufsize=1,
                )
        recvuntil("sim_output.log", match_string, timeout=60)
        return p


###


class TestLammpsIMDv3:

    @pytest.fixture()
    def command(self):
        return ["lmp"]

    @pytest.fixture()
    def input_file(self):
        # May cause issues, improve if needed
        absolute = Path("lammps/data/example_v3.in").resolve()
        return absolute

    @pytest.fixture()
    def match_string(self):
        return "Waiting for IMD connection on port 8888"

    @pytest.fixture()
    def simulation(self, tmp_path, command, match_string, input_file):
        old_cwd = Path.cwd()
        os.chdir(tmp_path)
        p = run_sim_and_wait(command, match_string, input_file=input_file)
        # can use path to check output files of simulation for validation
        yield tmp_path

        # Terminate the process
        p.terminate()
        try:
            p.wait(timeout=10)
        except subprocess.TimeoutExpired:
            p.kill()
            p.wait()
        finally:
            os.chdir(old_cwd)

    def test_compare_imd_to_dump(self, simulation):
        u = mda.Universe(
            Path(simulation / "topology_after_min.data"),
            "localhost:8888",
        )
        for ts in u.trajectory:
            pass
