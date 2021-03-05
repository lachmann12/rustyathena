use rand::Rng;
use frank::Ranking;

pub fn rank_rank_transform(arr: Vec<Vec<f32>>) -> Vec<Vec<i32>>{

    let mut rng = rand::thread_rng();

    let mut arr_rank = vec![vec![0f64; arr[0].len()]; arr.len()];
    for i in 0..arr.len() {
        arr_rank[i] = arr[i].rrank();
    }

    let mut arr = vec![vec![0f64; arr.len()]; arr[0].len()];
    for i in 0..arr_rank.len() {
        for j in 0..arr_rank[0].len() {
            arr[j][i] = arr_rank[i][j] + rng.gen::<f64>()/100000.0;
        }
    }

    let mut arr_rank = vec![vec![0i32; arr[0].len()]; arr.len()];
    for i in 0..arr.len() {
        let rr = arr[i].rrank();
        let mut rankround = vec![0i32; rr.len()];
        for j in 0..rr.len() {
            rankround[j] = rr[j] as i32;
        }
        arr_rank[i] = rankround;
    }
    return arr_rank;
}