"""DEPhT."""
from pathlib import Path
from flask import Blueprint, render_template, current_app, request

from depht.__main__ import run_depht

bp = Blueprint("depht", __name__)


def check_dir(dirpath):
    """Check if directory exists, otherwise make it."""
    if not dirpath.is_dir():
        dirpath.mkdir(parents=True)


@bp.route("/", methods=["GET"])
def depht():
    """Render depht."""
    length = current_app.config["MIN_SIZE"]
    cpu_cores = current_app.config["PHYSICAL_CORES"]
    att = current_app.config["ATT_SENSITIVITY"]
    local_models = current_app.config["LOCAL_MODELS"]

    '''
    # arguments to run depht
    # change this to a directory path??
    infile = request.args.get("input_file")
    outdir = Path("output")
    check_dir(outdir)

    model = request.args.get("model")
    runmode = request.args.get("mode")
    att_sens = request.args.get("att")
    no_draw = request.args.get("no_draw")
    dump_data = request.args.get("dump")
    products = request.args.get("products")
    length = request.args.get("length")
    cpus = request.args.get("cores")

    run_depht(infile, outdir, model=model, cpus=cpus, draw=no_draw,
              runmode=runmode, dump=dump_data, att_sens=att_sens,
              min_length=length)

    '''

    return render_template("depht/depht.html",
                           length=length,
                           cpu_cores=cpu_cores,
                           att=att,
                           local_models=local_models)


@bp.route("/results")
def depht_loading():
    """Render loading page for depht."""
    return render_template("loading_base.html")
