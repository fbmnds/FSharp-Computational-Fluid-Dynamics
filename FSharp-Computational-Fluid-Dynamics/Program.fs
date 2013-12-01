// Weitere Informationen zu F# unter "http://fsharp.net".
// Weitere Hilfe finden Sie im Projekt "F#-Lernprogramm".

open System
open System.Data
open System.IO
open Newtonsoft.Json
open Newtonsoft.Json.Linq
open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.LinearAlgebra.Double

open cfd_11

type Test_11_3 = { nx : int; 
                   dx : float; 
                   ny : int; 
                   dy : float; 
                   nt : int; 
                   dt : float; 
                   rho : float; 
                   nu : float; 
                   nit : int;
                   p0 : JArray;
                   u0 : JArray;
                   v0 : JArray;
                   u2 : JArray;
                   uu : JArray;
                   u_nt : JArray;
                   v_nt : JArray;
                   b0 : JArray;
                   b_nt : JArray;
                   p_nt : JArray;
                   p_py : JArray }

                   
let waitForKey() =
    printf "\n... enter key"
    System.Console.ReadKey() |> ignore


let getMatrix (x: JArray) nx ny : DenseMatrix = 
    let mutable m = DenseMatrix.create nx ny 0.
    // http://stackoverflow.com/questions/9976018/parsing-multidimensional-json-array-with-newtonsoft-json-net
    let mutable r = -1
    let mutable c = -1
    for i in x do
        r <- r + 1
        c <- -1
        for j in i do
            c <- c + 1
            m.At(r, c, (float j))
    m

let maxMatrix (x: DenseMatrix)  = 
    let mutable max = x.At(0,0)
    for i in [0..x.RowCount-1] do
        for j in [0..x.ColumnCount-1] do
            if x.At(i,j) > max then max <- x.At(i,j)
    max

let maxDiff (xj: JArray) (y: DenseMatrix) =
    let x = getMatrix xj y.RowCount y.ColumnCount
    let d = maxMatrix (((x - y) |> Matrix.map abs) :?> DenseMatrix)
    d


[<EntryPoint>]
let main argv = 

    let test_11_3 =
        __SOURCE_DIRECTORY__ + @"\test\test-11-3.json"
        |> File.ReadAllText

    let t = JsonConvert.DeserializeObject<Test_11_3> test_11_3
    // printfn "t.b0 : \n %A" (float (t.p0.First.First))


    let unt = getMatrix t.u_nt t.nx t.ny
    let vnt = getMatrix t.v_nt t.nx t.ny

    let bnt = getMatrix t.b_nt t.nx t.ny
    let bnt_ = buildUpB t.rho t.dt t.dx t.dy unt vnt

    //printfn "bnt_ : \n %A" bnt_
    
    printfn "max diff of t.b_nt, bnt_ : %.3f" (maxDiff t.b_nt bnt_)

    for i in [0..bnt_.RowCount-1] do
        for j in [0..bnt_.ColumnCount-1] do
            printf "%.5f " ((bnt_.At(i,j)) - (bnt.At(i,j)))
        printfn ""


    printfn "\n unt \n-----\n"
    for i in [0..bnt_.RowCount-1] do
        for j in [0..bnt_.ColumnCount-1] do
            printf "%.10f " (unt.At(i,j))
        printfn ""

    ////////////////////////
    /// verify pressPoission
    ////////////////////////

    let pnt = getMatrix t.p_nt t.nx t.ny
    let ppy = getMatrix t.p_py t.nx t.ny

    
    /// iterate 'nit' times: 
    let pP dx dy b p_ =
        let mutable p = p_
        for q in [0..t.nit-1] do
            p <- pressPoisson dx dy b p
        p

    let ppy_ = pP t.dx t.dy bnt pnt

    for i in [0..ppy_.RowCount-1] do
        for j in [0..ppy_.ColumnCount-1] do
            printf "%.7f " ((ppy_.At(i,j)) - (ppy.At(i,j)))
        printfn ""

    /////////////////////////
    /// => pressPoisson OK
    /////////////////////////

    ////////////////////////////////
    /// verify matrix multiplication
    ////////////////////////////////

    let u2 = getMatrix t.u2 t.nx t.ny
    let uu = getMatrix t.uu t.nx t.ny
    let u2_ = unt * unt

    /// Math.NET and Python differ:
    printfn "\nMath.NET and Python differ:"
    for i in [0..u2.RowCount-1] do
        for j in [0..u2.ColumnCount-1] do
            printf "%.7f " ((u2_.At(i,j)) - (u2.At(i,j)))
        printfn ""    

    /// Python is consistent: u*u = u**2
    printfn "\nPython is consistent: u*u = u**2:"
    for i in [0..u2.RowCount-1] do
        for j in [0..u2.ColumnCount-1] do
            printf "%.7f " ((uu.At(i,j)) - (u2.At(i,j)))
        printfn "" 

    waitForKey()
    0 // Exitcode aus ganzen Zahlen zurückgeben
