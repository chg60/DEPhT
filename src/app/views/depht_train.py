"""DEPhT Train."""
from flask import Blueprint, render_template, current_app, request
from depht_utils.pipelines import train_model

bp = Blueprint("depht_train", __name__)


@bp.route("/depht-train", methods=["GET"])
def depht_train():
    """Render depht."""
    window_size = current_app.config["WINDOW"]
    cpu_cores = current_app.config["LOGICAL_CORES"]

    # get arguments to train a new model

    model_name = request.args.get("name")
    phages = request.args.get("phage")
    bacteria = request.args.get("bact")
    prophage_csv = request.args.get("csv")
    window = request.args.get("window")
    cores = request.args.get("cpu")

    # train_model

    return render_template("depht_train/depht_train.html",
                           window=window_size,
                           cpu_cores=cpu_cores)
