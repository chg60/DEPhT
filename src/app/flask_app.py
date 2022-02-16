"""Flask app for DEPhT."""
from pathlib import Path

from flask import Flask

from app.views import depht, depht_train, documentation


PKG = "DEPhT"


def create_app(name=PKG, **kwargs):
    """
    Create the Flask app.

    :param name: name of the app
    :type name: str
    """
    app = Flask(name)
    app.config.from_object("config")

    path = Path(app.instance_path)

    if not path.is_dir():
        path.mkdir(parents=True)

    app.register_blueprint(depht.bp)
    app.register_blueprint(depht_train.bp)
    app.register_blueprint(documentation.bp)

    return app
