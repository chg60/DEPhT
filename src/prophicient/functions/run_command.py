from subprocess import Popen, DEVNULL, PIPE
import shlex


def run_command(command, verbose=False):
    """
    Runs the indicated command-line command, with output suppressed.

    `verbose`=False will discard any stdout and stderr output.
    `verbose`=True will return (stdout, stderr).

    :param command: command to run
    :type command: str
    :param verbose: dump stdout, stderr to console?
    :type verbose: bool
    :return: out, err
    """
    command = shlex.split(command)
    out, err = None, None
    if verbose:
        with Popen(args=command, stdout=PIPE, stderr=PIPE) as p:
            # Reading from stdout/stderr is blocking so don't need p.wait()
            out = p.stdout.read().decode("utf-8")
            err = p.stderr.read().decode("utf-8")
    else:
        with Popen(args=command, stdout=DEVNULL, stderr=DEVNULL) as p:
            # Not reading from stdout/stderr so we need to block using p.wait()
            p.wait()
    return out, err
