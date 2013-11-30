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



[<EntryPoint>]
let main argv = 
    printfn "%A" argv

    let test_11_3 =
        __SOURCE_DIRECTORY__ + @"\test\test-11-3.json"
        |> File.ReadAllText

    let t = JsonConvert.DeserializeObject<Test_11_3> test_11_3
    // printfn "t.b0 : \n %A" (float (t.p0.First.First))

    let u0 = getMatrix t.u0 t.nx t.ny
    let v0 = getMatrix t.v0 t.nx t.ny
    let b0 = getMatrix t.u0 t.nx t.ny
    let b0_ = buildUpB t.rho t.dt t.dx t.dy u0 v0

    printfn "b0 : \n %A" b0
    printfn "b0_ : \n %A" b0_
    
    let u_nt = getMatrix t.u_nt t.nx t.ny
    printfn "u_nt : \n %A" u_nt
    printfn "t.u_nt : \n %A" t.u_nt

    waitForKey()
    0 // Exitcode aus ganzen Zahlen zurückgeben
