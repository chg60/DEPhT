"""DEPhT."""
from pathlib import Path
from flask import Blueprint, render_template, current_app, request

from depht.__main__ import run_depht

bp = Blueprint("depht", __name__)


def check_dir(dirpath):
    """Check if directory exists, otherwise make it."""
    if not dirpath.is_dir():
        dirpath.mkdir(parents=True)


@bp.route("/")
def depht():
    """Render depht."""
    length = current_app.config["MIN_SIZE"]
    cpu_cores = current_app.config["LOGICAL_CORES"]
    att = current_app.config["ATT_SENSITIVITY"]
    local_models = current_app.config["LOCAL_MODELS"]

    return render_template("depht/depht.html",
                           length=length,
                           cpu_cores=cpu_cores,
                           att=att,
                           local_models=local_models)


@bp.route("/loading", methods=["GET", "POST"])
def depht_loading():
    """Render loading page for depht."""
    infile = [Path(request.args.get("input_file"))]
    outdir = Path("output")
    check_dir(outdir)

    args = []

    # this should be a str
    model = request.args.get("model")

    runmode = request.args.get("mode")
    att_sens = request.args.get("att")
    products = request.args.get("products")
    length = request.args.get("length")
    cpus = request.args.get("cores")

    # stash these as a Namespace

    # treat products as if None - replace with value later

    args = {"Input file": infile[0].name,
            "Model": model,
            "Runmode": runmode,
            "Att Sensitivity parameter": att_sens,
            "Min. phage homologs": products,
            "Min. length": length,
            "CPU cores": cpus}

    return render_template("depht/loading.html", args=args)


@bp.route("/results", methods=["GET", "POST"])
def results():
    """Display the results."""
    infile = [Path(request.args.get("input_file"))]
    outdir = Path("output")
    model = request.args.get("model")
    runmode = request.args.get("mode")
    att_sens = request.args.get("att")
    no_draw = request.args.get("no_draw")
    dump_data = request.args.get("dump")
    products = request.args.get("products")
    length = request.args.get("length")
    cpus = request.args.get("cores")

    print("starting depht")

    run_depht(infile, outdir, products, model=model, cpus=cpus, draw=no_draw,
              runmode=runmode, dump=dump_data, att_sens=att_sens,
              min_length=length)

    # write the output appropriately

    print("done running depht")

    return render_template("results.html")
