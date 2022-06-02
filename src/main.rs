// pub mod curve_v1;
pub mod curve_v1_1;



fn main(){
    let op: i128 = 150000000000;
    let ap: i128 = 40000000000;
    let of: i128 = 5500000000;

   

   let ask_amount = curve_v1_1::get_ask_amount_bisection(op, of, ap);

   println!("ask amount = {}, new ask pool = {}", ask_amount, ap - ask_amount);

   
   let get_offer: i128 = curve_v1_1::get_offer_amount_bisection(op, ask_amount, ap);
   println!("offer_amount = {}", get_offer);


    
}
