# Author: Magdalena Ondruskova

from typing import Tuple, Union

from pydantic import BaseModel


class RoutingCoordRequestBody(BaseModel):
    src_coord: Union[Tuple[float, float], None]
    dst_coord: Union[Tuple[float, float], None]
    from_time: str | None
    to_time: str | None
