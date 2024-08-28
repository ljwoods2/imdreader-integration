import imdreader
import MDAnalysis as mda
import pytest
import subprocess
import time
from pathlib import Path

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
    buffer = ""

    while time.time() < end_time:
        time.sleep(0.1)  # Small delay to avoid busy-waiting
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                buffer += line
                if target_line in line:
                    return line
    raise TimeoutError(
        f"Timeout after {timeout} seconds waiting for '{target_line}'"
    )


def run_sim_and_wait(command, match_string):
    with open("sim_output.log", "w") as f:
        p = subprocess.Popen(
            command,
            stdin=subprocess.PIPE,
            stdout=f,
            stderr=f,
            text=True,
            bufsize=1,
        )
        recvuntil("sim_output.log", match_string, timeout=10)
        return p


###


class TestLammpsIMDv3:

    @pytest.fixture()
    def command(self):
        # May cause issues, improve if needed
        absolute = Path("lammps/data/example_v3.in").resolve()
        return ["lmp", "<", absolute]

    @pytest.fixture()
    def match_string(self):
        return "Waiting for IMD connection on port 8888"

    @pytest.fixture()
    def simulation(self, tmp_path, command, match_string):
        with tmp_path.as_cwd():
            p = run_sim_and_wait(command, match_string)
            # can use path to check output files of simulation for validation
            yield tmp_path
            # Terminate the process
            p.terminate()
            try:
                p.wait(timeout=10)
            except subprocess.TimeoutExpired:
                p.kill()
                p.wait()

    def test_compare_imd_to_dump(self, simulation):
        u = mda.Universe(
            Path(simulation / "topology_after_min.data"),
            "localhost:8888",
        )
        for ts in u.trajectory:
            pass
