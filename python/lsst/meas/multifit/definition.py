from .multifitLib import definition_PositionComponent as PositionComponent
from .multifitLib import definition_RadiusComponent as RadiusComponent
from .multifitLib import definition_EllipticityComponent as EllipticityComponent
from .multifitLib import definition_Object as Object
from .multifitLib import definition_ObjectSet as ObjectSet
from .multifitLib import definition_Frame as Frame
from .multifitLib import definition_FrameSet as FrameSet
from .multifitLib import Definition

from . import multifitLib

multifitLib.Definition.Object = Object
multifitLib.Definition.Frame = Frame
multifitLib.Definition.ObjectSet = ObjectSet
multifitLib.Definition.FrameSet = FrameSet
multifitLib.Definition.PositionComponent = PositionComponent
multifitLib.Definition.RadiusComponent = RadiusComponent
multifitLib.Definition.EllipticityComponent = EllipticityComponent
