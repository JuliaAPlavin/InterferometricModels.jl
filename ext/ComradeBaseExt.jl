module ComradeBaseExt

import ComradeBase
import InterferometricModels
using InterferometricModels.StaticArrays

InterferometricModels.intensity(mod::ComradeBase.AbstractModel) =
    (x::InterferometricModels.XYType) -> ComradeBase.intensity_point(mod, (X=x[1], Y=x[2]))

InterferometricModels.visibility(mod::ComradeBase.AbstractModel) =
    (uv::InterferometricModels.UVType) -> ComradeBase.visibility_point(mod, (U=uv[1], V=uv[2]))

end
