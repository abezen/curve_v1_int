extern crate num;
use num::abs;


const A: i128 = 1000;
const PRECISION: i128 = 100;
const BETA: i128 = 10000;
const DIV: i128 = 10;
// const TEN_5: i128 = 100000;



pub fn get_sqrt(n: i128) -> i128 {
    let mut x: i128 = n;
    let mut y: i128 = 1;
    while x > y {
        x = (x + y)/2;
        y = n /x;
    }

    return x;
}

pub fn get_order(n: i128) -> i128 {
    let mut order = 1;
    let mut m: i128 = n /10;

    while m > 10 {
        m /= 10;
        order += 1;
    }

    return order;
}

pub fn get_step(n: i128) -> i128 {
    let mut num: i128 = 10;

    for _i in 1..n+1 {
        num *= 10;
    }

    return num;
}

pub fn get_initial_bisection_values_d(x0: i128, x1: i128) -> (i128, i128) {
    let d0: i128 = x0 + x1;

    let power: i128;

    if d0 < 100 {
        power = 0;
    } 
    else {
       let a: i128 = get_order(d0);
       power = a /3-1;
    }

    let d_step: i128 = get_step(power);
    let mut d_left: i128 = d0 - d_step;
    let mut d_right: i128 = d0 + d_step;
    let mut f_left: i128 = get_function_value_3(d_left, x0, x1);
    let mut f_right: i128 = get_function_value_3(d_right, x0, x1);

    
    while f_left > 0 && f_right > 0 {
        if f_left == f_right {
            d_left -= d_step;
            d_right += d_step;
            f_left = get_function_value_3(d_left, x0, x1);
            f_right = get_function_value_3(d_right, x0, x1);
        } else if f_left > f_right {
            d_right += d_step;
            f_right = get_function_value_3(d_right, x0, x1);
        } else {
            d_left -= d_step;
            f_left = get_function_value_3(d_left, x0, x1);
        }
    }

    while f_left < 0 && f_right < 0 {
        if f_left == f_right {
            d_left -= d_step;
            d_right += d_step;
            f_left = get_function_value_3(d_left, x0, x1);
            f_right = get_function_value_3(d_right, x0, x1);
        } else if f_left > f_right {
            d_left -= d_step;
            f_left = get_function_value_3(d_right, x0, x1);
        } else {
            d_right += d_step;
            f_right = get_function_value_3(d_left, x0, x1);
        }
    }

  //  println!("final d_left, d_right, f_left, f_right = {}, {}, {}, {}", d_left, d_right, f_left, f_right);
    
    return (d_left, d_right);
}

pub fn get_b(d: i128, x0: i128, x1: i128) -> i128 {
    let bk0: i128 = ( 4 * BETA * x0 * x1 + d * d / 2) / (d * d);
    let b: i128 = 1 + BETA - bk0;
    return b * b;  
}

pub fn get_function_value(d: i128, x0: i128, x1: i128) -> i128 {
    let b: i128 = get_b(d, x0, x1);
    let x01: i128 = x0 * x1;
    let t1: i128 = (4 * x01 * (x0 + x1 -d) + d /2)/d - b * (x01 - d * d);
    return t1 / b;
}

pub fn get_function_bisection_zero_d(x0: i128, x1: i128) -> i128 {
    let mut f_mid: i128;
    let (mut d_left, mut d_right) = get_initial_bisection_values_d(x0, x1);
    let mut d_mid: i128 = (d_left + d_right)/2;
    let mut f_left: i128 = get_function_value_3(d_left, x0, x1);
    let mut f_right: i128 = get_function_value_3(d_right, x0, x1);
    f_mid = get_function_value_3(d_mid, x0, x1);
    let mut d_mid_last: i128 = 0;

   // println!("d f_left init = {}, f_right init = {}", f_left, f_right);

   // while abs(f_mid) > PRECISION {
       while abs(d_mid - d_mid_last) > PRECISION {
           d_mid_last = d_mid.clone();

         //  println!("d mid last = {}, d left = {}, d right = {} ", d_mid_last, d_left, d_right);
         //  println!("f mid last = {}, f left = {}, f right = {} \n ", f_mid, f_left, f_right);
        if f_mid == 0  {
            return d_mid;
        }

        if f_left == 0 {
            return d_left;
        }

        if f_right == 0 {
            return d_right;
        }


        if f_left > 0 && f_mid > 0 && f_right < 0 && f_mid >= f_left {
            d_left = d_mid;
        }

    
        else if f_left < 0 && f_mid > 0 && f_right > 0 {
            d_right = d_mid;
        }

     
        else if f_left > 0 && f_mid < 0 && f_right < 0 {
            d_right = d_mid;
        }


        else if f_left < 0 && f_mid < 0 && f_right > 0 {
            d_left = d_mid;
        }

        
        f_left = get_function_value_3(d_left, x0, x1);
        f_right = get_function_value_3(d_right, x0, x1);
        d_mid = (d_left + d_right)/2;
        f_mid = get_function_value_3(d_mid, x0, x1);
    }
    println!("d = {}", d_mid);

    return d_mid;
}






pub fn get_function_value_3(da: i128, x0a: i128, x1a: i128) -> i128 {
    let pow: i128 = 10;
    let d: i128 =  da / pow;
    let x0: i128 = x0a / pow;
    let x1: i128 = x1a  / pow;

  //  let t1: i128 = 4 * A * (x0 + x1) + d - 4 * d * A;

    //if x0 == 0 || x1 == 0 {
    //    println!("x0 = {}, x1 = {}", x0, x1);
   // }
    //let t2: i128 = (d * d * d) / (4 * x0 * x1);

  //  if t1 == t2 {
  //      return 0;
  //  }

   // let f: i128 = t1 - t2;

   // if f == 0 {
      //  println!("x0 = {}, x1 = {}, d = {}", x0, x1, d);
   //     return 0;
  //  }
   
    return 4 * A * (x0 + x1) + d - 4 * d * A - (d * d * d) / (4 * x0 * x1);
    // return f;
    
}





pub fn get_ask_amount_bisection(op1: i128, of1: i128, ap1: i128) -> i128 {
    let op: i128 = op1 / DIV;
    let ap: i128 = ap1 / DIV;
    let of: i128 = of1 / DIV;

    println!("op = {}, ap = {}, of={} /div", op, ap, of);

    let d: i128 = get_function_bisection_zero_d(op, ap);

    println!("curve v1 d = {}", d);

    let new_ask_pool: i128 = get_function_bisection_zero_x(d, op+of);

    println!("new ask pool = {}", new_ask_pool);

    return (ap - new_ask_pool) * DIV;

}

pub fn get_offer_amount_bisection(op1: i128, aa1: i128, ap1: i128) -> i128 {
    let op: i128 = op1/DIV;
    let ap: i128 = ap1/DIV;
    let aa: i128 = aa1 / DIV;

    let d: i128 = get_function_bisection_zero_d(op, ap);

    let new_offer_pool: i128 = get_function_bisection_zero_x(d, ap - aa);
    return (new_offer_pool - op) * DIV;
}


pub fn get_function_bisection_zero_x(d: i128, x0: i128) -> i128 {
    
       let (mut x1_left, mut x1_right) = get_initial_bisection_values_x(d, x0);
       let mut x1_mid: i128 = (x1_left + x1_right)/2;
   
       let mut f_left: i128 = get_function_value_3(d, x0, x1_left);
       let mut f_right: i128 = get_function_value_3(d, x0, x1_right);
       let mut f_mid: i128 = get_function_value_3(d, x0, x1_mid);
        let mut x1_mid_last: i128 = 0;
   
       //while abs(f_mid) > PRECISION {
           while abs(x1_mid -x1_mid_last) > PRECISION {
               x1_mid_last = x1_mid;
          if f_mid == 0 {
            return x1_mid;
        }

        if f_left < 0 && f_mid < 0 && f_right > 0 {
            x1_left = x1_mid;
        }
        else if f_left < 0 && f_mid > 0 && f_right > 0 {
            x1_right = x1_mid;
        }
        else if f_left > 0 && f_mid > 0 && f_right < 0 {
            x1_left = x1_mid;
        }
        else if f_left > 0 && f_mid < 0 && f_right < 0 {
            x1_right = x1_mid;
        }

            x1_mid = (x1_left + x1_right)/2;
            f_mid = get_function_value_3(d, x0, x1_mid);
           f_left = get_function_value_3(d, x0, x1_left);
           f_right = get_function_value_3(d, x0, x1_right);
           
       }
   
       return x1_mid;
   }


   pub fn get_initial_bisection_values_x(d: i128, x0: i128 ) -> (i128, i128) {
   
   // let x1: i128 = d  / (2 *x0);

   let x1: i128 = d  - x0;
    let power: i128;

    if x1 < 100 {
        power = 0;
    } 
    else {
       let a: i128 = get_order(x1);
       power = a /3 - 1;
    }

    let x1_step: i128 = get_step(power);
    let mut x1_left: i128 = x1 - x1_step;
    let mut x1_right: i128 = x1 + x1_step;

    let mut f1_left: i128 = get_function_value_3(d, x0, x1_left);
    let mut f1_right: i128 = get_function_value_3(d, x0, x1_right);


    while f1_left == f1_right {
        x1_left -= x1_step;
        x1_right += x1_step;
        f1_left = get_function_value_3(d, x0, x1_left);
        f1_right = get_function_value_3(d, x0, x1_right);
    }

    if f1_left > 0 && f1_right > 0  && f1_left > f1_right {
        
        while  f1_right > 0 {
            x1_right += x1_step;
            f1_right = get_function_value_3(d, x0, x1_right);
        }
    }
    else if f1_left > 0 && f1_right > 0 && f1_left < f1_right  {
        
        while f1_left > 0 {
            x1_left -= x1_step;
            f1_left = get_function_value_3(d, x0, x1_left);
        }
    }
    else if f1_left < 0 && f1_right < 0 && f1_left > f1_right {
        
        while f1_left < 0 {
            x1_left -= x1_step;
            f1_left = get_function_value_3(d, x0, x1_left);
        }
    }
    else if f1_left < 0 && f1_right < 0 && f1_left < f1_right {
        
        while  f1_right < 0 {
            x1_right += x1_step;
            f1_right = get_function_value_3(d, x0, x1_right);
        }
    }

    return (x1_left, x1_right);
}






