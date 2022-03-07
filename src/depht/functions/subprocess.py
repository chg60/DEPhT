"""Wrapper around a subset of the subprocess and shlex built-in libraries
to make running subprocesses easy."""

import shlex
from subprocess import Popen, PIPE


def run_command(command):
    """Run the indicated command as a subprocess.

    If `verbose` is True, stdout and stderr will be returned and the
    code block that called this one can decide how to handle it.

    If the command invokes a program that is not globally executable,
    it will raise a FileNotFoundError with a message like "No such
    file or directory: 'xxx'".

    :param command: the command to run
    :type command: str
    :raises: FileNotFoundError
    :return: out, err
    """
    command = shlex.split(command)

    # Reading from stdout/stderr is blocking, no need to call sp.wait()
    with Popen(command, stdout=PIPE, stderr=PIPE, close_fds=True) as process:
        out = process.stdout.read().decode("utf-8")
        err = process.stderr.read().decode("utf-8")

    return out, err
