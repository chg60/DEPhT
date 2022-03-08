import pathlib
import shlex
from subprocess import DEVNULL, PIPE, Popen

# GLOBAL VARIABLES
# -----------------------------------------------------------------------------
MOL_TYPES = ["nucl", "prot"]


def create_blastdb(infasta, db_dir, db_name, mol_type="nucl", verbose=False,
                   parse_seqids=True, hash_index=False, gi_mask=False,
                   mask_data=None, mask_id=None, logfile=None,
                   tax_id=None, tax_id_map=None):
    if verbose:
        stdout = PIPE
    else:
        stdout = DEVNULL

    db_path = db_dir.joinpath(db_name)
    command = f"makeblastdb -in {infasta} -dbtype {mol_type} -out {db_path}"

    if parse_seqids:
        command = " ".join([command, "-parse_seqids"])
    if hash_index:
        command = " ".join([command, "-hash_index"])
    if gi_mask:
        command = " ".join([command, "-gi_mask"])
    if mask_data is not None:
        command = " ".join([command, "-mask_data", str(mask_data)])
    if mask_id is not None:
        command = " ".join([command, "-mask_id", str(mask_id)])
    if logfile is not None:
        command = " ".join([command, "-logile", str(logfile)])
    if tax_id is not None:
        command = " ".join([command,  "-taxid", str(tax_id)])
    if tax_id_map is not None:
        command = " ".join([command, "-taxid_map", str(tax_id_map)])

    split_command = shlex.split(command)
    with Popen(args=split_command, stdout=stdout, stderr=stdout) as process:
        out, errors = process.communicate()

    return db_path


def blastdbcmd(blast_db_path, entries, outfile_path, verbose=0):
    command = (f"blastdbcmd -db {blast_db_path} -out {outfile_path}")

    if isinstance(entries, list):
        command = " ".join([command, "-entry", ",".join([entries])])
    if isinstance(entries, str):
        command = " ".join([command, "-entry", entries])
    elif isinstance(entries, pathlib.Path):
        command = " ".join([command, "-entry_batch", str(entries)])

    command = shlex.split(command)
    with Popen(args=command, stdout=PIPE, stderr=PIPE) as process:
        out, errors = process.communicate()
        if verbose > 0 and out:
            print(out.decode("utf-8"))
        if verbose > 1 and errors:
            print(errors.decode("utf-8"))

