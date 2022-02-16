"""Documentation page."""
from flask import Blueprint, render_template

bp = Blueprint("documentation", __name__)


@bp.route("/documentation")
def documentation():
    """Render Home page."""
    return render_template("documentation.html")
