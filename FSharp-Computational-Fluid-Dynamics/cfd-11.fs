module cfd_11

open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra.Double


let SM_00 (z: DenseMatrix) = z.SubMatrix(0, (z.RowCount-2), 0, (z.ColumnCount-2)) :?> DenseMatrix
let SM_10 (z: DenseMatrix) = z.SubMatrix(1, (z.RowCount-2), 0, (z.ColumnCount-2)) :?> DenseMatrix
let SM_20 (z: DenseMatrix) = z.SubMatrix(2, (z.RowCount-2), 0, (z.ColumnCount-2)) :?> DenseMatrix
let SM_01 (z: DenseMatrix) = z.SubMatrix(0, (z.RowCount-2), 1, (z.ColumnCount-2)) :?> DenseMatrix
let SM_11 (z: DenseMatrix) = z.SubMatrix(1, (z.RowCount-2), 1, (z.ColumnCount-2)) :?> DenseMatrix
let SM_21 (z: DenseMatrix) = z.SubMatrix(2, (z.RowCount-2), 1, (z.ColumnCount-2)) :?> DenseMatrix
let SM_02 (z: DenseMatrix) = z.SubMatrix(0, (z.RowCount-2), 2, (z.ColumnCount-2)) :?> DenseMatrix
let SM_12 (z: DenseMatrix) = z.SubMatrix(1, (z.RowCount-2), 2, (z.ColumnCount-2)) :?> DenseMatrix
let SM_22 (z: DenseMatrix) = z.SubMatrix(2, (z.RowCount-2), 2, (z.ColumnCount-2)) :?> DenseMatrix


let buildUpB rho dt dx dy (u: DenseMatrix) (v: DenseMatrix) =
    let dimX2 = u.RowCount - 2
    let dimY2 = u.ColumnCount - 2
    let mutable bn = DenseMatrix.create u.RowCount u.ColumnCount 0.

    let x() : DenseMatrix = 
        let A = (SM_21 u) - (SM_01 u)
        let AA = A * A
        let B = (SM_12 v) - (SM_10 v)
        let BB = B * B
        let C = (SM_12 u) - (SM_10 u)
        let D = (SM_21 v) - (SM_01 v)
        let CD = C * D 
        let dxdt = rho/(2.*dx*dt)
        let dydt = rho/(2.*dy*dt)
        let dxdx = rho/(4.*dx*dx)
        let dxdy = rho/(2.*dx*dy)
        let dydy = rho/(4.*dy*dy)
        dxdt*A + dydt*B - dxdx*AA - dxdy*CD - dydy*BB
    bn.SetSubMatrix(1, dimX2, 1, dimY2, x())
    bn


let pressPoisson dx dy (b: DenseMatrix) (p: DenseMatrix) = 
    let mutable pn = p
    let dimX2 = p.RowCount - 2
    let dimY2 = p.ColumnCount - 2

    let x() : DenseMatrix =
        let A = (SM_21 p) + (SM_01 p)
        let B = (SM_12 p) + (SM_10 p)
        let C = (SM_11 b)
        let dx2 = dx*dx
        let dy2 = dy*dy
        let edxdy = 2. * (dx2 + dy2)
        let dxx = dx2 / edxdy
        let dyy = dy2 / edxdy
        dyy*A + dxx*B - dx2*dyy*C
    pn.SetSubMatrix(1, dimX2, 1, dimY2, x())

    pn.SetRow((dimX2+1), (pn.Row dimX2))
    pn.SetRow(0, (pn.Row 1))
    pn.SetColumn(0, (pn.Column 1))
    // redundant, if (because) the initial matrix p has already zeroed last column:
    pn.SetColumn((dimY2+1), (DenseVector.create pn.RowCount 0.)) 
    pn

