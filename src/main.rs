pub mod curve_v1;



fn main(){
    let op: i128 = 40000000000;
    let ap: i128 = 20000000000;
    let of: i128 = 1500000000;
    let d0: i128 = op + ap;

    let d = curve_v1::compute_d(op, ap, d0);

    let d_int: i128 = d as i128;

    println!("d = {}", d_int);
    

  // let ask_value = curve_v1::get_ask_amount_bisection(op, of, d, ap);

   let curve = curve_v1::curve_v1(op, ap, of);

   println!("The ask value = {}", curve);

   
   // let get_offer: i128 = curve_v1::compute_offer_amount_curve_v1(ap, op, ap);


    
}
