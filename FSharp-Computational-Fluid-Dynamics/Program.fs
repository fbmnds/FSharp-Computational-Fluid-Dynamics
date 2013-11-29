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


[<EntryPoint>]
let main argv = 
    printfn "%A" argv

    let test_11_3 =
        __SOURCE_DIRECTORY__ + @"\test\test-11-3.json"
        |> File.ReadAllText

    let t = JsonConvert.DeserializeObject<Test_11_3> test_11_3

    //let b0 = buildUpB t.rho t.dt t.dx t.dy t.u0 t.v0

    //printfn "b0 : \n %A" b0
    printfn "t.b0 : \n %A" (float (t.p0.First.First))
    
    // http://stackoverflow.com/questions/9976018/parsing-multidimensional-json-array-with-newtonsoft-json-net
    for i in (t.p0) do
        for j in i do
            printf " %.0f " (float j)
        printfn ""
    
    waitForKey()
    0 // Exitcode aus ganzen Zahlen zurückgeben
