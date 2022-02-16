"""About page."""
from flask import Blueprint, render_template

bp = Blueprint("about", __name__)


@bp.route("/about")
def about():
    """Render About page."""
    return render_template("about.html")
