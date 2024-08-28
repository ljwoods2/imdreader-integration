import imdreader
import MDAnalysis as mda
import pytest
import subprocess

def run_sim_and_wait(command, match_string):


    p = subprocess.Popen(
        command,
        stdin=subprocess.PIPE,
        stdout=f,
        stderr=f,
        text=True,
        bufsize=1,
    )
    return p
        


@pytest.fixture(scope="function")
def simulation(tmp_path):
    """Block until simulation is ready to accept connections.
    Yields the output directory for validation.
    Requires that the simulation engine is added to the PATH."""
    command = []
    match_string = "IMD: Will wait until I have a connection and IMD_GO orders."
    with tmp_path.as_cwd():
        p = run_sim_and_wait(command, match_string)
        yield tmp_path

        # Terminate the process
        p.terminate()
        try:
            p.wait(timeout=10)
        except subprocess.TimeoutExpired:
            p.kill()
            p.wait()


class TestLammpsIMDv3:
