#include <math.h>
#include <stdio.h>
#include <string>
#include <random>
#include <time.h>
#include <unordered_map>
#include "mex.h"
#include "matrix.h"
#include <chrono>

#define maxR 4096
#define maxN 524288
#define maxe 8192
#define maxK 512
#define accr 1e-12

void swap_double(double * a, int p, int q){
	double temp = a[p];
	a[p] = a[q];
	a[q] = temp;
}
void swap_int(int * a, int p, int q){
	int temp = a[p];
	a[p] = a[q];
	a[q] = temp;
}

void MYsort(double * a, double * b, int len, char type){
	// sort a according to rule type, permutate b according to a 
	if(len <= 1) return;
 	std::default_random_engine generator(time(0));
	std::uniform_int_distribution<int> distribution(0, len-1);
	int pos = distribution(generator), p = 0;
	double pivot = a[pos];
    swap_double(a, pos, len-1);
    swap_double(b, pos, len-1);
	if(type == 'd'){
		for(int i = 0; i < len-1; i++){
			if(a[i] > pivot){
				swap_double(a, i, p);
				swap_double(b, i, p);
			    p++;
			}
		}
		swap_double(a, p, len-1);
		swap_double(b, p, len-1);
	}
	if(type == 'a'){
		for(int i = 0; i < len-1; i++){
			if(a[i] < pivot){
				swap_double(a, i, p);
				swap_double(b, i, p);
			    p++;
			}
		}
		swap_double(a, p, len-1);
		swap_double(b, p, len-1);
	}
	MYsort(a, b, p, type);
	MYsort(a+p+1, b+p+1, len-p-1, type);
}

void projection_edge_weight(double * y, double * a, double * W, const double * para){
    y[0] = (W[0]*a[0] - W[1]*a[1])/(W[0] + W[1]);
    if(y[0] < 0) y[0] = fmax(-para[0], y[0]);
    else y[0] = fmin(para[0], y[0]);
    y[1] = - y[0];
}

void projection_para_edge_weight(double * y, double * a, double * W, const double * para, int len){
    int n = len/2, i, twoi;
    for(i = 0; i < n;i++){
        twoi = 2*i;
        y[twoi] = (W[twoi]*a[twoi] - W[twoi+1]*a[twoi+1])/(W[twoi] + W[twoi+1]);
        if(y[twoi] < 0) y[twoi] = fmax(-para[i], y[twoi]);
        else y[twoi] = fmin(para[i], y[twoi]);
        y[twoi+1] = - y[twoi];
    }
}


void projection_homo_weight(double * y, double * a, double * W, const double * para, int len){
	double * invW = new double[len];
    double suma = 0, suminvW = 0, lambda = 0, ypossum = 0, invWpossum = 0,
            ynegsum = 0, invWnegsum = 0, dlambda_thre = 0, tempy = 0;
    int *maskpos = new int [len];
    int *maskneg = new int [len];
    int maskposlen = 0, maskneglen = 0;
    int i;
    for(i = 0; i < len; i++){
        suma += a[i];
        invW[i] = 1/W[i];
        suminvW += invW[i];
        //mexPrintf("%f, %f", a[i], W[i]);
    }
    lambda = suma/suminvW;
    for(i = 0; i < len; i++){
        y[i] = -lambda * invW[i] + a[i];
        if(y[i] > accr){
            maskpos[maskposlen++] = i;
            ypossum += y[i];
            invWpossum += invW[i];
        }
        if(y[i] < -accr){
            maskneg[maskneglen++] = i;
            ynegsum += y[i];
            invWnegsum += invW[i];            
        }
    }
    while(ypossum > para[0] + accr && maskposlen >= 1){
        dlambda_thre = (ypossum - para[0])/invWpossum;
        for(i = maskposlen-1; i >=0; i--){
            //dlambda = min(y[maskpos[i]]*W[maskpos[i]],dlambda_thre);
            tempy = fmax(0, y[maskpos[i]]-dlambda_thre * invW[maskpos[i]]);
            ypossum -= y[maskpos[i]] - tempy;
            y[maskpos[i]] = tempy;
            if(y[maskpos[i]] < accr){
                invWpossum -= invW[maskpos[i]];
                swap_int(maskpos, i, maskposlen-1);
                maskposlen--;
            }
        }
    }
    while(ynegsum < -para[0] - accr && maskneglen >= 1){
        dlambda_thre = (-ynegsum - para[0])/invWnegsum;
        for(i = maskneglen-1; i >=0; i--){
            //dlambda = min(y[maskpos[i]]*W[maskpos[i]],dlambda_thre);
            tempy = fmin(0, y[maskneg[i]] + dlambda_thre * invW[maskneg[i]]);
            ynegsum += tempy - y[maskneg[i]];
            y[maskneg[i]] = tempy;
            if(y[maskneg[i]] > -accr){
                invWnegsum -= invW[maskneg[i]];
                swap_int(maskneg, i, maskneglen-1);
                maskneglen--;
            }
        }
    }
    delete[] invW;
    delete[] maskpos;
    delete[] maskneg;
}

void concave_card_subroutine_weight(double * x, double lambdamin, double lambdamax, double * a, double * W, const double * para, int i, int len){
	int j;
	double minval = INFINITY, tempsum = 0, lambda = 0, sumW =0;
    double * index = new double [len] , * val = new double [len], * invindex,
            * newa, * newW;
	int minindex;
//     printf("%d", len);
//     for(j = 0; j < len; j++){
//         printf("a[%d]:%f\n", j, a[j]);
//         printf("W[%d]:%f\n", j, W[j]);
//         printf("para[%d]:%f\n", j, para[j]);
//     }
    for(j = 0; j < len; j++){
        index[j] = j;
    }
	if(i%2 == 0){
		for(j = 0; j < len; j++){
            lambda -= para[j] + a[j];
            sumW += W[j];
        }
        lambda = lambda/sumW;
		for(j = 0; j < len; j++){
			val[j] = lambda*W[j] + a[j];
		}
        MYsort(val, index, len, 'a');
        for(j = 0; j < len; j++){
        	tempsum = tempsum + val[j] + para[j];
			if(tempsum < minval){
				minval = tempsum;
				minindex = j;
			}
        }
		if(minval > 0 || minindex == len - 1){
			for(j = 0; j < len; j++){
				x[j] = lambda;
			}
            delete[] index;
            delete[] val;
			return;
		}
	}else{
		lambda = (lambdamin + lambdamax)/2;
		for(j = 0; j < len; j++){
			val[j] = lambda*W[j] + a[j];
		}
        MYsort(val, index, len, 'a');
        for(j = 0; j < len; j++){
        	tempsum = tempsum + val[j] + para[j];
			if(tempsum < minval){
				minval = tempsum;
				minindex = j;
			}
        }       
		if(minval > 0){
			concave_card_subroutine_weight(x, lambdamin, lambda, a, W, para, i+1, len);
            delete[] index;
            delete[] val;
			return;
		}
		if(minindex == len-1){
			concave_card_subroutine_weight(x, lambda, lambdamax, a, W, para, i+1, len);
            delete[] index;
            delete[] val;
			return;
		}		
	}
    invindex = new double [len];
    newa = new double [len];
    newW = new double [len];
    for(j = 0; j < len; j++){
        invindex[(int)index[j]] = j;
        newa[j] = a[(int)index[j]];
        newW[j] = W[(int)index[j]];
    }
	concave_card_subroutine_weight(x, lambda, lambdamax, newa, newW, para, i+1, minindex+1);
	concave_card_subroutine_weight(x + minindex+1, lambdamin, lambda, newa + minindex + 1, newW + minindex + 1, para + minindex + 1, i+1, len - minindex - 1);
    double * newx = new double [len];
    for(j = 0; j < len; j++){
        newx[j] = x[(int)invindex[j]];
    }
    for(j = 0; j < len; j++){
        x[j] = newx[j];
    }    
    delete[] newx;
    delete[] index;
    delete[] invindex;
    delete[] newW;
    delete[] newa;
    delete[] val;

}

void projection_concave_card_weight(double * y, double * a, double * W, const double * para, int len){

    double * invW = new double [len];
    double * nega = new double [len];
    double lambdamax  = -INFINITY, lambdamin = INFINITY;
    int i;
	for(i = 0; i < len; i++){
		invW[i] = 1/W[i];
        nega[i] = -a[i];
        lambdamax = fmax(lambdamax, (-para[len-1]+a[i])*W[i]);
        lambdamin = fmin(lambdamin, (-para[0]+a[i])*W[i]);
	}    
    
	double * x = new double[len];
    // nega , invW may change
	concave_card_subroutine_weight(x, lambdamin, lambdamax, nega, invW, para, 0, len);
	for(i = 0; i < len; i++){
		y[i] = a[i] - x[i]/W[i];
	}
	delete[] invW;
	delete[] nega;
	delete[] x;
}

double gap_eff(int ** incidence_list, int * incidence_list_size,
        double ** parameter_list, int * parameter_list_size,
        char ** submodular_type, int * submodular_type_size,
        double * bias_vec, double * y_sum, double * x,
            int N, int R){
    double gap = 0, dv = 0, tempmax, tempmin;
    int i,j;
    for(i = 0; i < N; i++){
        gap -= y_sum[i]*x[i];
    }
    for(i = 0; i < R; i++){
        if(submodular_type[i][0] == 'h'){
            tempmax = -INFINITY;
            tempmin = INFINITY;
            for(j = 0; j < incidence_list_size[i]; j++){
                tempmax = fmax(tempmax, x[incidence_list[i][j]]);
                tempmin = fmin(tempmin, x[incidence_list[i][j]]);
            }
            gap += (tempmax - tempmin)*parameter_list[i][0];
        }
        if(submodular_type[i][0] == 'e'){
            gap += fabs(x[incidence_list[i][0]] - x[incidence_list[i][1]])*parameter_list[i][0];
        }
        if(submodular_type[i][0] == 'p'){
            for(j = 0; j < incidence_list_size[i]/2; j++){
                gap += fabs(x[incidence_list[i][2*j]] - x[incidence_list[i][2*j+1]])*parameter_list[i][j];
            }
        }
        if(submodular_type[i][0] == 'c'){
            double * tempx = new double [incidence_list_size[i]];
            double * index = new double [incidence_list_size[i]]();
            for(j = 0; j < incidence_list_size[i]; j++){
                tempx[j] = x[incidence_list[i][j]];
            }
            MYsort(tempx, index, incidence_list_size[i], 'd');
            for(j = 0; j < incidence_list_size[i]; j++){
                gap += tempx[j]*parameter_list[i][j];
            }
            delete[] tempx;
            delete[] index;
        }
    }
    
    return gap;
}

double gap_d_eff(int ** incidence_list, int * incidence_list_size,
        double ** parameter_list, int * parameter_list_size,
        char ** submodular_type, int * submodular_type_size,
        double * bias_vec, double * x,
            int N, int R, double * retarray){
    double * xsort = new double [N], * index = new double [N], * invindex = new double [N];
    double * level_cost = new double [N](), level_cost_min = 0;
    double pos_x_sum = 0, indexmax, indexmin;
    
    int i,j;
    double besti;
    for(i = 0; i < N; i++){
        xsort[i] = x[i];
        index[i] = i;
        if(x[i] > 0) pos_x_sum += x[i];
    }
    MYsort(xsort, index, N, 'd');
    for(i = 0; i < N; i++){
        invindex[(int)index[i]] = i;
    }
    for(i = 0; i < R; i++){
        if(submodular_type[i][0] == 'h'){
            indexmax = 0;
            indexmin = N-1;
            for(j = 0; j < incidence_list_size[i]; j++){
                indexmax = fmax(indexmax, invindex[incidence_list[i][j]]);
                indexmin = fmin(indexmin, invindex[incidence_list[i][j]]);
            }
            level_cost[(int)indexmin] += parameter_list[i][0];
            level_cost[(int)indexmax] -= parameter_list[i][0];
        }
        if(submodular_type[i][0] == 'e'){
            level_cost[(int)fmin(invindex[incidence_list[i][0]], invindex[incidence_list[i][1]])] += parameter_list[i][0];
            level_cost[(int)fmax(invindex[incidence_list[i][0]], invindex[incidence_list[i][1]])] -= parameter_list[i][0];
        }
        if(submodular_type[i][0] == 'p'){
            for(j = 0; j < incidence_list_size[i]/2; j++){
                level_cost[(int)fmin(invindex[incidence_list[i][2*j]], invindex[incidence_list[i][2*j+1]])] += parameter_list[i][j];
                level_cost[(int)fmax(invindex[incidence_list[i][2*j]], invindex[incidence_list[i][2*j+1]])] -= parameter_list[i][j];
            }
        }
        if(submodular_type[i][0] == 'c'){
            double * tempindex = new double [incidence_list_size[i]];
            double * nouseindex = new double [incidence_list_size[i]];
            for(j = 0; j < incidence_list_size[i]; j++){
                tempindex[j] = invindex[incidence_list[i][j]];
                nouseindex[j] = j;
            }
            MYsort(tempindex, nouseindex, incidence_list_size[i], 'a');
            for(j = 0; j < incidence_list_size[i]; j++){
                level_cost[(int)tempindex[j]] += parameter_list[i][j];
            }
            delete[] tempindex;
            delete[] nouseindex;
        }
    }
    for(i = 0; i < N; i++){
        level_cost[i] += bias_vec[(int)index[i]];
    }
    level_cost_min = fmin(level_cost[0], level_cost_min);
    for(i = 1; i < N; i++){
        level_cost[i] += level_cost[i-1];
        if(level_cost[i]< level_cost_min){
            level_cost_min = level_cost[i];
            besti = index[i];
        }           
    }

   
    delete[] index;
    delete[] invindex;
    delete[] xsort;
    delete[] level_cost;
//     mexPrintf("in %f\n", pos_x_sum);
//     mexPrintf("in %f\n", level_cost_min);   
    retarray[0] = level_cost_min;
    retarray[1] = pos_x_sum; 
    return besti;
}

void degree_comp(double ** degree_vec, int ** incidence_list, int * incidence_list_size, int * picked, int K){
    std::unordered_map<int, int> mymap;
    int i, j;
    for(i = 0; i < K; i++){
        for(j = 0; j < incidence_list_size[picked[i]]; j++){
            mymap[incidence_list[picked[i]][j]];
            mymap[incidence_list[picked[i]][j]]++;
        }
    }
    for(i = 0; i < K; i++){
        for(j = 0; j < incidence_list_size[picked[i]]; j++){
            degree_vec[i][j] = (double)mymap[incidence_list[picked[i]][j]];
        }
    }    
}
void complete_list(int * completelist, int * selected_list, int selected_list_size, int K, int R){
 	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, R-1);
    int count = selected_list_size, picked;
    int * flag = new int [R]();
    for(int i = 0; i < selected_list_size; i++){
        completelist[i] = selected_list[i];
        flag[selected_list[i]] = 1;
    }
    while(count < K){
        picked = distribution(generator);
        if(flag[picked]!=1){
            completelist[count] = picked;
            flag[picked] = 1;
            count++;
        }
    }
    delete[] flag;
}

void main_func(double * y_sum, double * record, double * record_discrete,
        int ** incidence_list, int * incidence_list_size,
        double ** parameter_list, int * parameter_list_size,
        char ** submodular_type, int * submodular_type_size,
        double * bias_vec, int N, int R, int K, unsigned int T,
        unsigned int record_dis, int G_part,
        int ** greedy_part_list, int * greedy_part_size, double gaptol, double * timevec,
        double * Xmat, double * sweepi){
    double ** y = new double * [R];
    int * completelist = new int [K];
    double ** newy = new double * [K];
    double ** a = new double * [K];
    double ** degree_vec = new double * [K];
    double * x = new double [N]();
    int i, j, t, selected;
    double tempeta, eta;
    
    // New things added
    auto start_time = std::chrono::high_resolution_clock::now();
    auto current_time = std::chrono::high_resolution_clock::now();
    double time2;
    double level_cost_min;
    double pos_x_sum;
    double * retarray = new double [2]();
    double besti;
    int place, ii;
    // End new declarations
    
 	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, G_part-1);
    //int picked[maxK]; 
    double para = 1, gap;
    for(i = 0; i < R; i++){
        y[i] = new double [incidence_list_size[i]]();
    }
    for(i = 0; i < K; i++){
        newy[i] = new double [maxe]; 
        a[i] = new double [maxe]; 
        degree_vec[i] = new double [maxe]; 
    }
    for(t = 0; t < T; t++){
        selected = distribution(generator);
        // to complete each part to be with size K
        complete_list(completelist, greedy_part_list[selected], greedy_part_size[selected], K, R);
        degree_comp(degree_vec, incidence_list, incidence_list_size, completelist, K);
        for(j = 0; j < K; j++){
            for(i = 0; i < incidence_list_size[completelist[j]]; i++){
                a[j][i] = y[completelist[j]][i] -
                        (y_sum[incidence_list[completelist[j]][i]]+ bias_vec[incidence_list[completelist[j]][i]])/degree_vec[j][i];
                //mexPrintf("%f\n", a[j][i]);
                //mexPrintf("a:%f, %d,%d\n", picked[j], i, a[picked[j]][i]);
                
            }
            if(submodular_type[completelist[j]][0] == 'h'){
                projection_homo_weight(newy[j], a[j], degree_vec[j], parameter_list[completelist[j]], incidence_list_size[completelist[j]]);
            }
            if(submodular_type[completelist[j]][0] == 'e'){
                projection_edge_weight(newy[j], a[j], degree_vec[j], parameter_list[completelist[j]]);
            }
            if(submodular_type[completelist[j]][0] == 'p'){
                projection_para_edge_weight(newy[j], a[j], degree_vec[j], parameter_list[completelist[j]], incidence_list_size[completelist[j]]);
            }
            if(submodular_type[completelist[j]][0] == 'c'){
                projection_concave_card_weight(newy[j], a[j], degree_vec[j], parameter_list[completelist[j]], incidence_list_size[completelist[j]]);
            }
        }
        for(j = 0; j < K; j++){
            for(i = 0; i < incidence_list_size[completelist[j]]; i++){
                y_sum[incidence_list[completelist[j]][i]] += newy[j][i] - y[completelist[j]][i];
                y[completelist[j]][i] = newy[j][i];
            }
        }
        if((t+1)%record_dis == 0){
            for(i = 0; i < N; i++){
                x[i] = - y_sum[i] - bias_vec[i];
            }
            // Focus on discrete gap
            current_time = std::chrono::high_resolution_clock::now();
            time2 = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
            mexPrintf("Time: %f\n", time2);
            
            besti = gap_d_eff(incidence_list, incidence_list_size,
                parameter_list, parameter_list_size,
                submodular_type, submodular_type_size,
                bias_vec, x, N, R, retarray);
            level_cost_min = retarray[0];
            pos_x_sum = retarray[1];           
            gap = level_cost_min+pos_x_sum;
            mexPrintf("%f\n", pos_x_sum);
            mexPrintf("%f\n", level_cost_min);
//             mexPrintf("gap = %f\n", gap);
            place = (t+1)/record_dis-1;
            
            // Save the time, the level_cost_min, and pos_x_sum
            // to later re-create gap if desired, and have level_costs
           
            record_discrete[place] = level_cost_min;
            record[place] = pos_x_sum;
            timevec[place] = time2;
            sweepi[place] = besti;
            for(ii = 0; ii < N; ii++){
                Xmat[place*N+ii] = x[ii];
            }
             if(gap < gaptol){
               	return;
            }
            //End changes            
        }
    }
    delete[] y;
    delete[] newy;
    delete[] a;
    delete[] degree_vec;
    delete[] x;
    delete[] completelist;
    delete[] retarray;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	if (nlhs != 6 || nrhs != 12) {
    mexWarnMsgTxt("Check Parameters");
    return;
    }
    // read incidence_list
    int R = *(mxGetPr(prhs[5]));
    const mxArray * incidence_list_org = prhs[0];
    mxArray * incidence_Element;
    //const int * R_org = mxGetDimensions(incidence_list_org);
    int * * incidence_list = new int * [R];
    int * incidence_list_size = new int [R];
    double * templist;
    int j, k;
    for(j = 0; j < R;j++){
       incidence_Element = mxGetCell(incidence_list_org, j);
       incidence_list_size[j] = (int)mxGetN(incidence_Element);
       incidence_list[j] = new int [incidence_list_size[j]];
       templist = mxGetPr(incidence_Element);
       for(k = 0; k < incidence_list_size[j]; k++){
           incidence_list[j][k] = templist[k] - 1;
       }
    }
    // read parameter_list
    const mxArray * parameter_list_org = prhs[1];
    mxArray * parameter_Element;
    double * * parameter_list = new double * [R];
    int * parameter_list_size = new int [R];
    for(j = 0; j < R;j++){
       parameter_Element = mxGetCell(parameter_list_org, j);
       parameter_list_size[j] = (int)mxGetN(parameter_Element);
       parameter_list[j] = new double[parameter_list_size[j]];
       templist = mxGetPr(parameter_Element);
       for(k = 0; k < parameter_list_size[j]; k++){
           parameter_list[j][k] = templist[k];
           //mexPrintf("%d, %d, %f\n", j, k, parameter_list[j][k]);
       }
    }
    // read submodular_type
    const mxArray * submodular_type_org = prhs[2];
    mxArray * submodular_type_Element;
    char ** submodular_type = new char * [R];
    int * submodular_type_size = new int [R];
    char * temp_type;
    for(j = 0; j < R;j++){
       submodular_type_Element = mxGetCell(submodular_type_org, j);
       temp_type = (char *)mxGetPr(submodular_type_Element);
       submodular_type_size[j] = (int)mxGetN(submodular_type_Element);
       submodular_type[j] = new char[submodular_type_size[j]];
       for(k=0; k < submodular_type_size[j]; k++){
           submodular_type[j][k] = (char)temp_type[2*k];
           //mexPrintf("%d, %d, %c\n", j, k, temp_type[2*k]);
       }
    }
    // read greedy_part
    int G_part = *(mxGetPr(prhs[9]));
    const mxArray * greedy_part = prhs[10];
    mxArray * greedy_part_Element;
    int ** greedy_part_list = new int * [G_part];
    int * greedy_part_size = new int [G_part];
    double * greedy_part_templist;
    for(j = 0; j < G_part;j++){
       greedy_part_Element = mxGetCell(greedy_part, j);
       greedy_part_size[j] = (int)mxGetN(greedy_part_Element);
       greedy_part_list[j] = new int [greedy_part_size[j]];
       greedy_part_templist = mxGetPr(greedy_part_Element);
       for(k = 0; k < greedy_part_size[j]; k++){
           greedy_part_list[j][k] = greedy_part_templist[k] - 1;
           //mexPrintf("%d, %d, %d", greedy_part_list[j][k], j, k);
       }
    }
        
    double * bias_vec = mxGetPr(prhs[3]);
    int N = *(mxGetPr(prhs[4]));
    int K = *(mxGetPr(prhs[6]));
    unsigned int T = *(mxGetPr(prhs[7]));
    unsigned int record_dis = *(mxGetPr(prhs[8]));
    double gaptol = *(mxGetPr(prhs[11]));

    int record_num = T/record_dis;
    
    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    double * y_sum = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(record_num, 1, mxREAL);
    double * record = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(record_num, 1, mxREAL);
    double * record_discrete = mxGetPr(plhs[2]);
    //New left hand side output
    plhs[3] = mxCreateDoubleMatrix(record_num, 1, mxREAL);
    double * timevec = mxGetPr(plhs[3]);
    
    //New output mat
    plhs[4] = mxCreateDoubleMatrix(1,record_num*N, mxREAL);
    double * Xmat = mxGetPr(plhs[4]);
//     mexPrintf("length(X) = %d",record_num*N);

    //Best i for sweep cut
    plhs[5] = mxCreateDoubleMatrix(record_num, 1, mxREAL);
    double * sweepi = mxGetPr(plhs[5]);
    
    main_func(y_sum, record, record_discrete,
            incidence_list, incidence_list_size,
            parameter_list, parameter_list_size,
            submodular_type, submodular_type_size,
            bias_vec, N, R, K, T, record_dis, G_part,
            greedy_part_list, greedy_part_size, gaptol,timevec,
            Xmat, sweepi);
    delete[] incidence_list;
    delete[] incidence_list_size;
    delete[] parameter_list;
    delete[] parameter_list_size;
    delete[] submodular_type;
    delete[] submodular_type_size;
    delete[] greedy_part_list;
    delete[] greedy_part_size;
}







