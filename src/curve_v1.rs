extern crate num;

//use num_integer::Roots;
// use cosmwasm_std::{Decimal256, Uint256};

const A: i128 = 85;

const LEFT: i128 = 0;
const RIGHT: i128 = 100000000000;
const PRECISION: i128 = 10;



 use num::abs;
/* ----------------------------------------------------------- */
/// A trait for things that can be approximately equal.
pub trait Epsilon {
    type RHS;
    type Precision;

    /// Return true if self and `other` differ no more than by a given amount.
    fn close(&self, other: Self::RHS, precision: Self::Precision) -> bool;

    /// Return true if self is close to zero.
    fn near_zero(&self, precision: Self::Precision) -> bool;
}

impl Epsilon for i128 {
    type RHS = i128;
    type Precision = i128;

    fn close(&self, other: i128, precision: i128) -> bool {
        abs(other - *self) < abs(precision)
    }

    fn near_zero(&self, precision: i128) -> bool {
        abs(*self) < abs(precision)
    }
}


/* ---------- Newton's method (single root) ---------- */

/// Configuration structure for the Newton's method (one root version).
#[derive(Debug, Clone, Copy)]
pub struct OneRootNewtonCfg {
    /// The real root, if any, is most likely to be within this distance from
    /// the reported root, but this is not guaranteed.
    pub precision: i128,
    /// A limit on the number of iterations to perform. Pass `None` if you
    /// don't want a limit.
    pub max_iters: Option<u32>
}

pub fn newton_one(config: OneRootNewtonCfg,
                           first_approx: i128,
                           b: i128,
                           q: i128,
                           s: i128
                          )
    -> Option<i128>
   // where T: Float + Epsilon<RHS=T, Precision=T>
          
{
    let mut left = LEFT;
    let mut right = RIGHT;
    let mut left_val = get_d_target(left, b, q, s);
    let mut right_val = get_d_target(right, b, q, s);
    let mut root = first_approx;
    let mut prev_root = None;
    let mut iter = 0;
    while prev_root.map_or(true, |old| !root.close(old, config.precision))
        && config.max_iters.map_or(true, |max| iter < max) {
            iter += 1;
            if let Some(next) = next_newton_iter(config.precision,
                                                root,
                                                 b,
                                                 q,
                                                 s,
                                                 left,
                                                 right
                                                 ) {
                prev_root = Some(root);
                root = next;
            } else if let Some(fallback_root) 
                = linear_fallback(left, right, left_val, right_val) {
                    prev_root = Some(root);
                    root = fallback_root;
            } else {
                return None
            }
            let val_at_root = get_d_target(root, b, q, s);
            if left_val * val_at_root <= 0 {
                right = root;
                right_val = val_at_root;
            } else {
                left = root;
                left_val = val_at_root;
            }
    }
    Some(root)
}

fn next_newton_iter(prec: i128,
                             old: i128,
                             b: i128,
                             q: i128,
                             s: i128,
                             left: i128,
                             right: i128
                             )
    -> Option<i128>     
{
    let d = get_d_deriv(old, b, q );
    if d.near_zero(prec) {
        return None
    }
   
   let res = get_d_function(old, b, q, s);
    if res < left {
        None
    } else if res > right {
        None
    } else {
        Some(res)
    }
}

fn linear_fallback(x1: i128 , x2: i128, y1: i128, y2: i128) -> Option<i128>
{
    let res = ((y2 - y1) * x1 - (x2 - x1) * y1) / (y2 - y1);
    if res < x1 {
        None
    } else if res > x2 {
        None
    } else {
        Some(res)
    }
}

pub fn get_d_function(x: i128, b: i128, q: i128, s: i128) -> i128 {
    return (2 * x * x * x + 4 * A *s * q + (3 * x * x + b * q)/2) / (3 * x * x + b * q);
}

pub fn get_d_target(x: i128, b: i128, q: i128, s: i128) -> i128 {
        return (x * x * x + x * b * q - 4 * A * s * q) / q;
}

pub fn get_d_deriv(x: i128, b: i128, q: i128) -> i128 {
    return (3 * x * x + b * q) / q;
}

/* ---------------------------------
    IMPLEMENTATION OF THE CURVE V1 ALGORITHM ACCORDING TO 
    THE PAPER BY MICHAEL EGOROV ON 10 NOVEMBER 2019
*/

pub fn compute_d(op: i128, ap: i128, d0: i128) -> i128
{
    let sum: i128 = op + ap;
    let prod: i128 =  op * ap;
    let a4: i128 = 4 * A;
    let prod4: i128 = 4 * prod;
    let a4_1: i128 = a4 - 1;
  
    let _cfg = OneRootNewtonCfg {
        precision: PRECISION,
        max_iters: None
    };

   
    let sol = newton_one(_cfg, d0, a4_1, prod4, sum);

    let d: i128;

    match sol {
        Some(ss) => d = ss,
        None => panic!(),
    };

    return d;
}

pub fn get_x_target(x: i128, x1:i128, d: i128) -> i128 {

    println!("get_x_target A, x, x1, d = {}, {}, {}, {}", A, x, x1, d);
    return (16*A*d*x1*x + d*d*d-16*A*x1*x*(x1+x)-4*x1*x*d)/(4*x1*x);
}

pub fn get_x_deriv(x: i128, x1: i128, d: i128) -> i128 {
    return (d*d*d+16*A*x1*x*x)*(-1)/(4*x1*x*x);
}


pub fn get_x_function(x: i128, x1: i128, d: i128) -> i128 {
    let num1: i128 = -2 *d*d*d-16*A*d*x1*x+4*d*x1*x+16*A*x1*x1*x;
  
    let num2: i128 = -d*d*d-16*A*x1*x*x;
   
  

    let n1: i128 = num1 / num2;
    let n2: i128 = num1 % num2;

   

    let res: i128 = n1 * x + (n2 * x) / num2;

    return res;
}




pub fn curve_v1(_offer_pool: i128, _ask_pool: i128, _offer: i128)  -> i128
{
    let op = _offer_pool as i128;
    let ap = _ask_pool as i128;
    let of = _offer as i128;
    //let d0 = 1000.0;
    let d0: i128 = op + ap;
    let d = compute_d(op, ap, d0);

    //println!("d version 2 = {0}", d);
    let ask_amnt: i128 = get_ask_amount_bisection(op, of, d, ap);
  //   let ask_f = get_ask_amount(op, of, d, d0);
   // let ask_amnt: i128 = ask_f as i128;
    return _ask_pool - ask_amnt;
}






pub fn compute_offer_amount_curve_v1(ask_pool: i128, offer_pool: i128, ask_amount: i128)  -> i128
{   
    let op = offer_pool as i128;
    let ap = ask_pool as i128;
    let am = ask_amount as i128;
  
     let d0: i128 = op + ap;

    let d = compute_d(op, ap, d0);
  
   let offer_amnt: i128 = get_ask_amount_bisection(ap, am, d, op);
  
   
    return offer_amnt - offer_pool;
}


pub fn get_function_value(d: i128, x1: i128, x: i128) -> i128 {
    let num: i128 = 16 * A * d * x1 * x + d* d*d -16 * A * (x1 +x) * x1 * x - 4 * d * x1 * x;
    let denom = 4 * x1 *x;

    return num / denom;
    
    
}

pub fn get_ask_amount_bisection(op: i128, of: i128, d:i128, ap: i128) -> i128 {
    
    let x1: i128 = op + of;

    let mut t1: i128 = ap - of;
    let mut t2: i128 = ap + of;

    let mut t_mid: i128 = 0;
    let mut y_mid;
    
    
    
    let mut y1 = get_function_value(d, x1, t1);
    let mut y2 = get_function_value(d, x1, t2);

    while y1 * y2 > 0 {
        t1 -= of;
        t2 += of;

        y1 = get_function_value(d, x1, t1);
        y2 = get_function_value(d, x1, t2);

    }

    
    while abs(y1 -y2) > 1000 {

        y1 = get_function_value(d, x1, t1);
        y2 = get_function_value(d, x1, t2);
        t_mid = (t1 + t2)/2;
        y_mid = get_function_value(d, x1, t_mid);

        if y1 * y_mid <= 0 {
            t2 = t_mid;
        } 
        else {
            t1 = t_mid;
        }
    }
    

    return t_mid;
}