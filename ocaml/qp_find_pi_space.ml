open Core
open Qputils
open Qptypes

let run ?(sym="None") ?(apply="no") fc ezfio_filename =
  Ezfio.set_file ezfio_filename ;

  let aos = 
    match Input.Ao_basis.read () with
    | Some aos -> aos
    | None -> failwith "Unable to read AOs"
  and mos = 
    match Input.Mo_basis.read () with
    | Some mos -> mos
    | None -> failwith "Unable to read MOs"
  in
  let rec find_power symmetry accu = function
    | -1 -> accu
    | i  -> let new_accu = 
              if aos.Input.Ao_basis.ao_power.(i) = symmetry then (i::accu) 
              else accu
           in
           find_power symmetry new_accu (i-1)
  and n =
    (AO_number.to_int aos.Input.Ao_basis.ao_num) 
  and m =
    (MO_number.to_int mos.Input.Mo_basis.mo_tot_num) 
  and one  = Positive_int.of_int 1
  and zero =  Positive_int.of_int 0
  in

  (* Indices of the px, py and pz AOs *)
  let x_indices =
    find_power Symmetry.Xyz.{x=one  ; y=zero ; z=zero} [] (n-1)
  and y_indices =
    find_power Symmetry.Xyz.{x=zero ; y=one  ; z=zero} [] (n-1)
  and z_indices =
    find_power Symmetry.Xyz.{x=zero ; y=zero ; z=one } [] (n-1)
  in

  (* Compute the relative weight of each MO on the px, py, pz spaces *)
  let compute_weight mo_i list_aos =
    let num =
       List.fold_left ~f:(fun s i -> s +. (MO_coef.to_float @@ mos.Input.Mo_basis.mo_coef.(mo_i).(i)) ** 2.) ~init:0. list_aos
    and denom =
       Array.fold ~f:(fun s x -> s +. (MO_coef.to_float x) ** 2.) ~init:0. mos.Input.Mo_basis.mo_coef.(mo_i) 
    in
    num /. denom
  in
  let result = 
    Array.init ~f:(fun mo_i -> 
      (mo_i+1,
      compute_weight mo_i x_indices,
      compute_weight mo_i y_indices,
      compute_weight mo_i z_indices) ) m
    |> Array.to_list
  in

  let pi, sigma =
    let rec aux test_xyz (accu_pi: int list) (accu_sigma: int list) = function
    | [] -> (List.rev accu_pi, List.rev accu_sigma)
    | (mo_i,x,y,z)::rest ->
        if test_xyz (x,y,z) then
           aux test_xyz (mo_i::accu_pi) accu_sigma rest
         else
           aux test_xyz accu_pi (mo_i::accu_sigma) rest
    in
    match sym with
    | "x" | "X" -> aux (fun (x,y,z) -> (x>y && x>z)) [] [] result
    | "y" | "Y" -> aux (fun (x,y,z) -> (y>x && y>z)) [] [] result
    | "z" | "Z" -> aux (fun (x,y,z) -> (z>x && z>y)) [] [] result
    | _         -> ([],[])
  in

  begin
    match sym with
    | "x" | "X"  | "y" | "Y"  | "z" | "Z" -> 
        if apply = "cas" then
          begin
            let ne = Ezfio.get_electrons_elec_alpha_num () in
            let command = 
              "qp_set_mo_class -act \"[" ^
              (String.concat ~sep:"," @@ List.map ~f:string_of_int pi) ^
              "]\" -core \"[" ^
              (String.concat ~sep:"," @@ List.map ~f:string_of_int @@ List.filter ~f:(fun mo_i -> mo_i <= ne) sigma) ^
              "]\" -del \"[" ^
              (String.concat ~sep:"," @@ List.map ~f:string_of_int @@ List.filter ~f:(fun mo_i -> mo_i > ne) sigma) ^
              "]\" " ^ ezfio_filename
            in
            print_endline command;
            if Sys.command command <> 0 then
              failwith "Command failed"
          end
        else if apply = "cassd" then
          begin
            let ne = Ezfio.get_electrons_elec_alpha_num () in
            let nfc = 
              if fc then
                let ic =
                  Unix.open_process_in ("qp_set_frozen_core.py -q "^ezfio_filename)
                in
                let result = int_of_string @@ input_line ic in
                close_in ic; result
              else 0
            in
            let command = 
              "qp_set_mo_class" ^
              ( if nfc > 0 then
                Printf.sprintf " -core \"[1-%d]\"" nfc else "") ^
              " -inact \"[" ^
              (String.concat ~sep:"," @@ List.map ~f:string_of_int @@ List.filter ~f:(fun mo_i -> mo_i > nfc && mo_i <= ne) sigma) ^ "]\"" ^
              " -act \"[" ^ (String.concat ~sep:"," @@ List.map ~f:string_of_int pi) ^ "]\"" ^
              " -virt \"[" ^
              (String.concat ~sep:"," @@ List.map ~f:string_of_int @@ List.filter ~f:(fun mo_i -> mo_i > ne) sigma) ^ "]\" " ^ ezfio_filename
            in
            print_endline command;
            if Sys.command command <> 0 then
              failwith "Command failed"
          end
        else
          begin
            Printf.printf "Pi: [";
            List.iter ~f:(fun mo_i -> Printf.printf "%d," mo_i) pi;
            Printf.printf "\b]\n\nSigma: [";
            List.iter ~f:(fun mo_i -> Printf.printf "%d," mo_i) sigma;
            Printf.printf "\b]\n"
          end
    | _ -> List.iter ~f:(fun (mo_i,x,y,z) -> Printf.printf "%d: (%f,%f,%f)\n" mo_i x y z) result
  end

  



let spec =
  let open Command.Spec in
  empty
  +> flag "sym" (optional string)
    ~doc:"{x,y,z} Axis perpendicular to the plane"
  +> flag "apply" (optional string)
    ~doc:"[cas|cassd] Set the pi space as the active space, all other MOs frozen. If [cassd] is chosen, set inactive and virtual mos."
  +> flag "fc" no_arg
    ~doc:"freeze core MOs with -apply cassd"
  +> anon ("ezfio_filename" %: string)



let command =
    Command.basic_spec
    ~summary: "Quantum Package command"
    ~readme:(fun () ->
     "Find all the pi molecular orbitals to create a pi space.
     ")
    spec
    (fun sym apply fc ezfio_filename () -> run ?sym ?apply fc ezfio_filename)


let () =
    Command.run command

