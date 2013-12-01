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
        let AA = A.PointwiseMultiply(A) :?> DenseMatrix
        let B = (SM_12 v) - (SM_10 v)
        let BB = B.PointwiseMultiply(B) :?> DenseMatrix
        let C = (SM_12 u) - (SM_10 u)
        let D = (SM_21 v) - (SM_01 v)
        let CD = C.PointwiseMultiply(D) :?> DenseMatrix 
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


let cavityFlow rho nu nit nt dt dx dy nx ny p_ u_ v_ =
    let dtdx = dt/dx
    let dtdy = dt/dy
    let dtrhodx = dt/(2.*rho*dx)
    let dtrhody = dt/(2.*rho*dy)
    let nudtdx2 = (nu*dt)/(dx*dx)
    let nudtdy2 = (nu*dt)/(dy*dy)
    let mutable p = p_
    let mutable u = u_
    let mutable v = v_
    let mutable b = buildUpB rho dt dx dy u v
    for i1 in [0..nt-1] do
        if i1 > 0 then b <- buildUpB rho dt dx dy u v
        for i2 in [0..nit-1] do
            p <- pressPoisson dx dy b p
        let Bu = (SM_10 u)
        let Du = (SM_01 u)
        let Eu = (SM_11 u)
        let Fu = (SM_21 u)
        let Hu = (SM_12 u)
        let Bv = (SM_10 v)
        let Dv = (SM_01 v)
        let Ev = (SM_11 v)
        let Fv = (SM_21 v)
        let Hv = (SM_12 v)
        let Bp = (SM_10 p)
        let Dp = (SM_01 p)
        let Fp = (SM_21 p)
        let Hp = (SM_12 p)
        let EuEuDu = Eu.PointwiseMultiply(Eu - Du) :?> DenseMatrix
        let EvEuBu = Ev.PointwiseMultiply(Eu - Bu) :?> DenseMatrix
        let FpDp = Fp - Dp
        let Fu2EuDu = Fu - 2.*Eu + Du
        let Hu2EuBu = Hu - 2.*Eu + Bu
        let HpBp = Hp - Bp
        let EuEvDv = Eu.PointwiseMultiply(Ev - Dv) :?> DenseMatrix
        let EvEvBv = Ev.PointwiseMultiply(Ev - Bv) :?> DenseMatrix
        let Fv2EvDv = Fv - 2.*Ev + Dv
        let Hv2EvBv = Hv - 2.*Ev + Bv
        let ux = Eu - dtdx*EuEuDu - dtdy*EvEuBu - dtrhodx*FpDp + nudtdx2*Fu2EuDu + nudtdy2*Hu2EuBu
        let vx = Ev - dtdx*EuEvDv - dtdy*EvEvBv - dtrhody*HpBp + nudtdx2*Fv2EvDv + nudtdy2*Hv2EvBv
        u.SetSubMatrix(1, (u.RowCount-2), 1, (u.ColumnCount-2), ux)
        v.SetSubMatrix(1, (v.RowCount-2), 1, (v.ColumnCount-2), vx)
        u.SetRow(0, (DenseVector.create u.RowCount 0.))
        u.SetColumn(0, (DenseVector.create u.ColumnCount 0.))
        u.SetColumn(u.ColumnCount-1, (DenseVector.create u.ColumnCount 1.)) // ## in last line overwritten below
        v.SetRow(0, (DenseVector.create v.RowCount 0.))
        v.SetRow(v.RowCount-1, (DenseVector.create v.RowCount 0.))
        v.SetColumn(0, (DenseVector.create v.ColumnCount 0.))
        v.SetColumn(v.ColumnCount-1, (DenseVector.create v.ColumnCount 0.))
        u.SetRow(u.RowCount-1, (DenseVector.create u.RowCount 0.))
    (u,v,p,b)