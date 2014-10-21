open Core.Std;;

let input_data = "
* Positive_float : float  
  assert (x >= 0.) ; 

* Strictly_positive_float : float  
  assert (x > 0.) ; 

* Negative_float : float  
  assert (x <= 0.) ; 

* Strictly_negative_float : float  
  assert (x < 0.) ; 

* Positive_int : int  
  assert (x >= 0) ; 

* Strictly_positive_int : int  
  assert (x > 0) ; 

* Negative_int : int  
  assert (x <= 0) ; 

* Det_coef : float
  assert (x >= -1.) ; 
  assert (x <=  1.) ; 

* Normalized_float : float
  assert (x <= 1.) ; 
  assert (x >= 0.) ; 

* Strictly_negative_int : int  
  assert (x < 0) ; 

* Non_empty_string : string
  assert (x <> \"\") ;

* Atom_number : int  
  assert (x > 0) ; 
  if (x > 1000) then
    warning \"More than 1000 atoms\";
  if (Ezfio.has_nuclei_nucl_num ()) then
    assert (x <= (Ezfio.get_nuclei_nucl_num ()));

* MO_number : int  
  assert (x > 0) ; 
  if (x > 1000) then
    warning \"More than 1000 MOs\";
  if (Ezfio.has_mo_basis_mo_tot_num ()) then
    assert (x <= (Ezfio.get_mo_basis_mo_tot_num ()));

* AO_number : int  
  assert (x > 0) ; 
  if (x > 1000) then
    warning \"More than 1000 AOs\";
  if (Ezfio.has_ao_basis_ao_num ()) then
    assert (x <= (Ezfio.get_ao_basis_ao_num ()));

* Nucl_number : int  
  assert (x > 0) ; 
  if (x > 1000) then
    warning \"More than 1000 Atoms\";
  if (Ezfio.has_nuclei_nucl_num ()) then
    assert (x <= (Ezfio.get_nuclei_nucl_num ()));

* N_int_number : int 
  assert (x > 0) ; 
  if (x > 100) then
    warning \"N_int > 100\";
  if (Ezfio.has_determinants_n_int ()) then
    assert (x = (Ezfio.get_determinants_n_int ()));

* Det_number : int 
  assert (x > 0) ; 
  if (x > 100000000) then
    warning \"More than 100 million determinants\";
  if (Ezfio.has_determinants_det_num ()) then
    assert (x <= (Ezfio.get_determinants_det_num ()));

* States_number : int 
  assert (x > 0) ; 
  if (x > 100) then
    warning \"More than 100 states\";
  if (Ezfio.has_determinants_n_states_diag ()) then
    assert (x <= (Ezfio.get_determinants_n_states_diag ()))
  else if (Ezfio.has_determinants_n_states ()) then
    assert (x <= (Ezfio.get_determinants_n_states ()));

* Bit_kind_size : int  
  begin match x with
  | 8 | 16 | 32 | 64 -> ()
  | _ -> raise (Failure \"Bit_kind_size should be (8|16|32|64).\")
  end;

* Bit_kind : int  
  begin match x with
  | 1 | 2 | 4 | 8 -> ()
  | _ -> raise (Failure \"Bit_kind should be (1|2|4|8).\")
  end;

* MO_coef : float

* AO_coef : float

"
;;

let untouched = "
module Determinant : sig
  type t
  val to_int64_array : t -> int64 array
  val of_int64_array : int64 array -> t
  val to_string : t -> string
end = struct
  type t = int64 array
  let to_int64_array x = x
  let of_int64_array x = 
    if (Ezfio.has_determinants_n_int ()) then
      begin
        let n_int = Ezfio.get_determinants_n_int () in
        assert ((Array.length x) = n_int*2) 
      end
      ; x
  let to_string x = Array.to_list x
    |> List.map ~f:Int64.to_string
    |> String.concat ~sep:\", \"
end

"

let template = format_of_string "
module %s : sig
  type t
  val to_%s : t -> %s
  val of_%s : %s -> t
  val to_string : t -> string
end = struct
  type t = %s
  let to_%s x = x
  let of_%s x = ( %s x )
  let to_string x = %s.to_string x
end

"
;;

let parse_input input=
  print_string "open Core.Std;;\nlet warning = print_string;;\n" ;
  let rec parse result = function
    | [] -> result
    | ( "" , ""   )::tail -> parse result tail
    | ( t  , text )::tail -> 
        let  name , typ  = String.lsplit2_exn ~on:':' t
        in
        let typ  = String.strip typ
        and name = String.strip name in
        let typ_cap = String.capitalize typ in
        let newstring = Printf.sprintf template name typ typ typ typ typ typ typ 
          ( String.strip text ) typ_cap
        in
        List.rev (parse (newstring::result) tail )
  in
     String.split ~on:'*' input
  |> List.map ~f:(String.lsplit2_exn ~on:'\n') 
  |> parse []
  |> String.concat 
  |> print_string
;;

let () = 
  parse_input input_data ;
  print_endline untouched
;;

