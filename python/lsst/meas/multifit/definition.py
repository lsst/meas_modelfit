from .multifitLib import definition_PositionElement as PositionElement
from .multifitLib import definition_RadiusElement as RadiusElement
from .multifitLib import definition_EllipticityElement as EllipticityElement
from .multifitLib import definition_ObjectComponent as ObjectComponent
from .multifitLib import definition_ObjectComponentSet as ObjectComponentSet
from .multifitLib import definition_Frame as Frame
from .multifitLib import definition_FluxGroup as FluxGroup
from .multifitLib import definition_FrameSet as FrameSet
from .multifitLib import Definition

from . import multifitLib

multifitLib.Definition.ObjectComponent = ObjectComponent
multifitLib.Definition.FluxGroup = FluxGroup
multifitLib.Definition.Frame = Frame
multifitLib.Definition.ObjectComponentSet = ObjectComponentSet
multifitLib.Definition.FrameSet = FrameSet
multifitLib.Definition.PositionElement = PositionElement
multifitLib.Definition.RadiusElement = RadiusElement
multifitLib.Definition.EllipticityElement = EllipticityElement

PositionElement.Bounds = multifitLib.detail_CircleConstraint
RadiusElement.Bounds = multifitLib.detail_MinMaxConstraint
EllipticityElement.Bounds = multifitLib.detail_CircleConstraint
