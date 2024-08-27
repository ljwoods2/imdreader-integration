import imdreader
import MDAnalysis as mda
import pytest
import subprocess

def run_sim_and_wait(command, match_string):
    """Block until simulation is ready to accept connections.
    Yields the output directory for validation.
    Requires that the simulation engine is added to the PATH."""
    with tmpdir.as_cwd():
        with open("gmx_output.log", "w") as f:
            p = subprocess.Popen(
                command,
                stdin=subprocess.PIPE,
                stdout=f,
                stderr=f,
                text=True,
                bufsize=1,
            )
            try:
                yield port
            finally:
                # Terminate the process
                p.terminate()
                try:
                    p.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    p.kill()
                    p.wait()

                # Ensure all file descriptors are closed
                f.close()


@pytest.fixture(scope="function")
def simulation():
    command = []
    match_string = "IMD: Will wait until I have a connection and IMD_GO orders."
    run_sim_and_wait(command, match_string)