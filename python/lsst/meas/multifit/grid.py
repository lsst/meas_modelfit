from .multifitLib import grid_PositionElement as PositionElement
from .multifitLib import grid_RadiusElement as RadiusElement
from .multifitLib import grid_EllipticityElement as EllipticityElement
from .multifitLib import grid_ObjectComponent as ObjectComponent
from .multifitLib import grid_SourceComponent as SourceComponent
from .multifitLib import grid_Frame as Frame
from .multifitLib import grid_FluxGroup as FluxGroup
from .multifitLib import Grid

from . import multifitLib

multifitLib.Grid.ObjectComponent = ObjectComponent
multifitLib.Grid.SourceComponent = SourceComponent
multifitLib.Grid.Frame = Frame
multifitLib.Grid.FluxGroup = FluxGroup
multifitLib.Grid.PositionElement = PositionElement
multifitLib.Grid.RadiusElement = RadiusElement
multifitLib.Grid.EllipticityElement = EllipticityElement
multifitLib.Grid.PositionArray = multifitLib.grid_PositionArray
multifitLib.Grid.RadiusArray = multifitLib.grid_RadiusArray
multifitLib.Grid.EllipticityArray = multifitLib.grid_EllipticityArray
multifitLib.Grid.ObjectComponentArray = multifitLib.grid_ObjectComponentArray
multifitLib.Grid.SourceComponentArray = multifitLib.grid_SourceComponentArray
multifitLib.Grid.FluxGroupArray = multifitLib.grid_FluxGroupArray
multifitLib.Grid.FrameArray = multifitLib.grid_FrameArray

FluxGroup.ComponentArray = multifitLib.grid_FluxGroup_ComponentArray
ObjectComponent.SourceArray = multifitLib.grid_ObjectComponent_SourceArray
