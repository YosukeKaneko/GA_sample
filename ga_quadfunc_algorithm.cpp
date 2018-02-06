#include <stdio.h>
#include <string.h>
#include <iostream>
#include <cmath>
#include <bitset>
#include <random>
#include <string>

using namespace std;

// 2次関数クラス
class quad_function
{
public:
    double a;
    double b;
    double c;

    quad_function(double a, double b, double c)
    {
        this->a = a;
        this->b = b;
        this->c = c;
    }

    float calc (float x)
    {
        return a*pow(x, 2) + b * x + c;
    }

};

class binary_to_decimal
{
public:
    int GENE_LENGTH;
    double min, max;

    binary_to_decimal(int GENE_LENGTH, double min, double max)
    {
        this->GENE_LENGTH = GENE_LENGTH;
        this->min = min;
        this->max = max;
    }

    // TODO：変換する桁数を可変にする 
    double convert(string binary)
    {
        unsigned long decimal = bitset<16>(binary).to_ulong();
        return min + (decimal / (pow(2, GENE_LENGTH) - 1) ) * (max - min);
    }
};

// void initialize_gene_pool(string gene_pool[100][16], int GENE_LENGTH, int GENE_POOL_LENGTH)
void initialize_gene_pool(string gene_pool[100], int GENE_LENGTH, int GENE_POOL_LENGTH)
{
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_int_distribution<> rand100(0, pow(2, GENE_LENGTH)); // [0, 2^n] 範囲の一様乱数

    for(int i = 0; i < GENE_POOL_LENGTH; i++)
    {
        // TODO：変換する桁数を可変にする 
        string random_init_num = bitset<16>(rand100(mt)).to_string();
        gene_pool[i] = random_init_num;
        // cout << random_init_num << endl;
        // for(int j = 0; j < GENE_LENGTH; j++)
        // {   
        //     gene_pool[i][j] = random_init_num[j];
        // }
    }
}

double set_random_score(double max)
{
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_real_distribution<double> rand100(0.0, max);// [0, max] 範囲の一様乱数
    return rand100(mt);
}
int set_random_int(int max)
{
    random_device rnd;     // 非決定的な乱数生成器を生成
    mt19937 mt(rnd());     //  メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    uniform_int_distribution<int> rand100(0, max);// [0, max] 範囲の一様乱数
    return rand100(mt);
}

// 先頭から順に2つの遺伝子ペアを作り、ある確率（probability変数で設定）で交叉を実行する
// todo gene_poolの配列サイズを動的に設定できるようにしたい
void simple_crossover(string gene_pool[100],int GENE_LENGTH, int GENE_POOL_LENGTH, double probability)
{
    for(int i = 0; i < GENE_POOL_LENGTH; i += 2)
    {
        // probabilityの確率で交叉を実行
        if(set_random_score(100.0) > probability) { continue; }
        // 交叉する位置をランダムに決定
        int cross_over_point = set_random_int(GENE_LENGTH - 1);
        for(int j = 0; j < GENE_LENGTH; j++)
        {
            if(j <= cross_over_point)
            {
                gene_pool[i][j] = gene_pool[i + 1][j];
            }
            else
            {
                gene_pool[i + 1][j] = gene_pool[i][j];
            }
        }
    }
}

// 各個体の全遺伝子が、ある確率（probability変数で設定）で突然変異を起こす
// todo gene_poolの配列サイズを動的に設定できるようにしたい
void simple_mutation(string gene_pool[100], int GENE_LENGTH, int GENE_POOL_LENGTH, double probability)
{
    for(int i = 0; i < GENE_POOL_LENGTH; i++)
    {
        for(int j = 0; j < GENE_LENGTH; j++)
        {
            if(set_random_score(100.0) > probability){ continue; }
            // gene_pool[i][j] = to_string(set_random_int(1));
        }
    }
}

double find_answer(quad_function quadfunc, double MIN, double MAX)
{
    double max_value = 0;
    double temp_x = MIN;
    int STEPS = 100;
    double step = (MAX - MIN) / STEPS;
    for(int i = 0; i < STEPS; i++)
    {
        max_value = (max_value < quadfunc.calc(temp_x))? quadfunc.calc(temp_x) : max_value;
        temp_x += step;
    }
    return max_value;
}

int main()
{
    // 遺伝子のbit数と、プールサイズの定義
    static const int GENE_LENGTH = 16;
    static const int GENE_POOL_LENGTH = 100;
    // 評価関すのパラメーターを定義：f(x) = Ax^2 + Bx + C
    static const double A = -1;
    static const double B = 1;
    static const double C = 0;
    // 探索する地域を設定
    static const double MIN = 0.0;
    static const double MAX = 1.0;
    // 最大の世代数を定義
    static const int MAX_GENERATIONS = 100000;

    // その他必要な一時保存用の変数を宣言
    // double eval_val[GENE_POOL_LENGTH];
    double arrow_score;
    double temp_score;
    double total_score;
    int selected_gene_index[GENE_POOL_LENGTH];
    string best_gene;
    double highest_score;
    double highest_input_decimal;
    string old_gene_pool[GENE_POOL_LENGTH];


    // 2次関数インスタンスをnewする
    quad_function quadfunc(A, B, C);
    // cout << quadfunc.calc(2);
    
    // 2進数->10進数変換インスタンスをnewする
    binary_to_decimal btod(GENE_LENGTH, MIN, MAX);
    // cout << btod.convert("1110111111111111");

    // 遺伝子の初期化    
    string gene_pool[GENE_POOL_LENGTH] = {string(GENE_LENGTH, '0')};
    initialize_gene_pool(gene_pool, GENE_LENGTH, GENE_POOL_LENGTH);

    // 世代交代ループ
    for(int gen = 0; gen < 1; gen++)
    {
        // 遺伝子ごとの評価値の計算
        total_score = 0;

        for(int i = 0; i < GENE_POOL_LENGTH; i++)
        {
            best_gene = ( highest_score <= quadfunc.calc(btod.convert(gene_pool[i])) ) ? gene_pool[i] : best_gene;
 
            total_score += quadfunc.calc(btod.convert(gene_pool[i]));
            // cout << eval_val[i] << endl;
        }
        // 遺伝子の選択
        for(int i = 0; i < GENE_POOL_LENGTH; i++)
        {
            temp_score = 0;
            arrow_score = set_random_score(total_score);

            for(int j = 0; j < GENE_POOL_LENGTH; j++)
            {
                temp_score += quadfunc.calc(btod.convert(gene_pool[j]));
                if(temp_score >=  arrow_score){
                    selected_gene_index[i] = j;
                    break;
                }
            }
        }

        // 選択された遺伝子indexの表示
        // for(int i = 0; i < 100 ; i++)
        // {
        //     cout << selected_gene_index[i] << endl;
        // }

        // 新しい遺伝子プールを作るために遺伝子プールをコピー（関数で書き換えたい）
        // memcpy(old_gene_pool, gene_pool, GENE_POOL_LENGTH);
        for(int i = 0; i < GENE_POOL_LENGTH ; i++)
        {
            old_gene_pool[i] = gene_pool[i];
        }
        // 遺伝子プールの選択
        for(int i = 0; i < GENE_POOL_LENGTH; i++)
        {
            gene_pool[i] = old_gene_pool[selected_gene_index[i]];
        }

        // 交叉の実行
        simple_crossover(gene_pool, GENE_LENGTH, GENE_POOL_LENGTH, 80.0);

        // 突然変異の実行
        // simple_mutation(gene_pool, GENE_LENGTH, GENE_POOL_LENGTH, 20.0); 

    }

    // cout << "最高評価値：" << highest_score << endl;
    cout << "最高評価値の遺伝子："  << best_gene << endl;
    cout << "最高評価値の遺伝子の10進数変換："  << btod.convert(best_gene) << endl;
    cout << "最高評価値：" << quadfunc.calc(btod.convert(best_gene)) << endl;
    cout << endl;
    cout << "正解の評価値：" << find_answer(quadfunc, MIN, MAX) << endl;

    return 0;
}