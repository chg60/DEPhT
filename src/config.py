"""Config file."""

from depht.__main__ import (MIN_CDS_FEATURES, ATT_SENSITIVITY,
                            MIN_SIZE, MIN_PRODUCTS_NORMAL, MIN_PRODUCTS_STRICT,
                            PHYSICAL_CORES, LOCAL_MODELS)
from depht_utils.pipelines.train_model import WINDOW, LOGICAL_CORES

MIN_CDS_FEATURES = MIN_CDS_FEATURES
ATT_SENSITIVITY = ATT_SENSITIVITY

MIN_SIZE = MIN_SIZE
MIN_PRODUCTS_NORMAL = MIN_PRODUCTS_NORMAL
MIN_PRODUCTS_STRICT = MIN_PRODUCTS_STRICT

WINDOW = WINDOW
LOGICAL_CORES = LOGICAL_CORES

PHYSICAL_CORES = PHYSICAL_CORES

# check that .zip files etc are not added
LOCAL_MODELS = [
    model for model in LOCAL_MODELS if not len(model.split(".")) > 1]
