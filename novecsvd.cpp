#include<bits/stdc++.h>
// #include "SVD.h"

using namespace std;
const int N=29;
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX(x,y) ((x)>(y)?(x):(y))
int SVD(float a[][N], int m, int n,float *w, float v[][N]);
bool sortinrev(const pair < float,vector<float> > &a,  
               const pair < float,vector<float> > &b) 
{ 
       return (a.first > b.first); 
}
int main()
{
    time_t start, end; 
    start = clock();
    ios_base::sync_with_stdio(false); 

    int rw =24998;
    int row = 29;
    int col = 29;
    // vector< vector<float> > r(row, vector<float> (col,0));
    float r[row][N] ={0.0};
    // r[0][0] = 1;
    // r[0][1] = 4;
    // r[0][2] = 0;
    // r[1][0] = 1;
    // r[1][1] = 1; 
    // r[1][2] = 1;
    // r[2][0] = 3;
    // r[2][1] = 0;
    // r[2][2] = 4;
    float data[rw][col];
    std::ifstream file("intrusion.csv");
    //float a = 0.0;
    std::string lin;
    std::getline(file, lin);
    //cout<<lin;
    for(int rows = 0; rows < rw; ++rows)
    {
        std::string line;
        std::getline(file, line);
        if ( !file.good() )
            break;
        std::stringstream iss(line);
        for (int cols = 0; cols < col; ++cols)
        {
            std::string val;
            std::getline(iss, val, ',');
            if ( !iss.good() )
                break;

            std::stringstream convertor(val);
            convertor >> data[rows][cols];
            // if(rows==0)
            // cout<<data[rows][cols]<<" ";
        }
        //cout<<endl;
    }
    // vector<float> mean(col,0);
    // vector<float> std(col,0);
    float mean[col] = {0};
    float std[col] ={0};

    for(int i = 0;i<col;i++)
    {
        float t = 0;
        for(int j = 0;j<rw;j++)
        {
            t = t + data[j][i];
        }
        t= t/rw;
        mean[i] = t;
    }

    for(int i = 0;i<col;i++)
    {
        float t = 0;
        for(int j = 0;j<rw;j++)
        {
            float t2 = data[j][i]- mean[i];
            t = t + t2*t2;
        }
        t= t/rw;
        std[i] = sqrt(t);
    }
    // for(int i = 0;i<col;i++)
    // {

    //     cout<<mean[i]<<"mean"<<endl;

    // }

    // for(int i = 0;i<col;i++)
    // {

    //     cout<<std[i]<<"std"<<endl;

    // }

    for(int i = 0;i<rw;i++)
    {
        for(int j = 0;j<col;j++)
        {
            data[i][j] = (data[i][j]- mean[j])/std[j];
        }
    }

    

    for(int i =0;i<col;i++)
    {
        for(int j=0;j<col;j++)
        {
            for(int k=0;k<rw;k++)
            {
                r[i][j] = r[i][j] + data[k][i]*data[k][j];
            }
        }
    }
    // cout<<"{{{{{{";
    // for(int i = 0;i<col;i++)
    // {
    //     cout<<data[0][i]<<endl;
    // }

    // for(int i = 0;i<col;i++)
    // {
    //     cout<<r[0][i]<<endl;
    // }

    // vector<float> w(row);
	// vector< vector<float> > v(row, vector<float> (col,0));
    float w[row];
    float v[row][N] ={0};
   
    SVD(r, row, col, w, v);
    // for(int i =0;i<1;i++)
    // {
    //     cout<<"[";
    //     for(int j = 0;j<col;j++)
    //     {
    //         cout<<v[i][j]<<" , ";
    //     }
    //     cout<<"]";
    //     cout<<endl;
    // }

    //******
    //sort(w.begin(),w.end());
    float eigen_sum  = 0;
    for(int i =0;i<row;i++)
    {    
        eigen_sum  += fabs(w[i]);
        // cout<<"*****"<<w[i]<<"\n ";   
    }
    // cout<<endl;

    float thres= 0.9;

    vector< vector<float> > v1(row, vector<float> (col,0));

    for(int i =0;i<row;i++)
    {
        for(int j=0;j<col;j++)
        {
            v1[i][j] = v[i][j];
        }
    }

    vector< pair < float,vector<float> > > vc;

    for(int i =0;i<row;i++)
    {
        vc.push_back( make_pair(w[i],v1[i]) );

    }
    sort(vc.begin(), vc.end(),sortinrev);

    float sum = 0;
    int dim = 0;
    for(int i =0;i<row;i++)
    {
        sum += fabs(vc[i].first);
        dim++;
        float temp = sum/eigen_sum;
        if(temp>0.9)
        {
            break;
        }

    }
    // cout<<endl;
    // cout<<"DIMENSION  "<<eigen_sum<<"  "<<dim<<endl;

    // cout<<endl;
    // for(int i =0;i<col;i++)
    // cout<<vc[0].second[i]<<"yes"<<endl;
    ///cout<<endl;
    
    //cout<<vc[0].second[0]<<"yes"<<endl;


    vector< vector<float> > red_data(rw, vector<float> (dim,0));


    for(int i =0;i<rw;i++)
    {
        for(int j=0;j<dim;j++)
        {
            for(int k=0;k<col;k++)
            {
                red_data[i][j] = red_data[i][j] + data[i][k]*(vc[j].second[k]);
            }
        }
    }

    // for(int i =0;i<2;i++)
    // {
    //     for(int j=0;j<dim;j++)
    //     {
    //         cout<<red_data[i][j];
    //         cout<<endl;
    //     }
    // }

      

    ofstream outfile;
    outfile.open ("reduced_data.csv");
    for (int i = 0;i<dim-1;i++)
    {
        outfile<<"col"<<i<<",";
    }
    outfile<<"col"<<dim-1<<endl;

    for(int i =0;i<rw;i++)
    {
        for(int j=0;j<dim-1;j++)
        {
            outfile<<red_data[i][j]<<",";
        }
        outfile<<red_data[i][dim-1]<<endl;
    }

    outfile.close();
    end = clock();  
  
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    cout << "Time taken by program is : " << fixed  
         << time_taken << setprecision(5); 
    cout << " sec " << endl;
    //myfile << "Writing this to a file.\n";
    // myfile.close();
    // for(int i =0;i<row;i++)
    // {
    //     cout<<"[";
    //     for(int j = 0;j<col;j++)
    //     {
    //         cout<<r[i][j]<<" , ";
    //     }
    //     cout<<"]";
    //     cout<<endl;
    // }
}

// int SVD(vector< vector<float> > &m,vector< vector<float> > &v,vector<float>&w,int row,int col)
// {
// 	//vector<float> w(row);
// 	//vector<vector<float>> v(row, vector<int> (col,0));
// 	SVD(m, row, col, w, v);
// 	//m.SortCols(w);//****implement this *****
// 	return 0;
// }

///////////////////////////////////////////////////////

static double PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}

int SVD(float a[][N], int m, int n,float *w, float v[][N])
{
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < n) 
    {
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }
  
    rv1 = (double *)malloc((unsigned int) n*sizeof(double));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++) 
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m) 
        {
            for (k = i; k < m; k++) 
                scale += fabs((double)a[k][i]);
            if (scale) 
            {
                for (k = i; k < m; k++) 
                {
                    a[k][i] = (float)((double)a[k][i]/scale);
                    s += ((double)a[k][i] * (double)a[k][i]);
                }
                f = (double)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (float)(f - g);
                if (i != n - 1) 
                {
                    for (j = l; j < n; j++) 
                    {
                        for (s = 0.0, k = i; k < m; k++) 
                            s += ((double)a[k][i] * (double)a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a[k][j] += (float)(f * (double)a[k][i]);
                    }
                }
                for (k = i; k < m; k++) 
                    a[k][i] = (float)((double)a[k][i]*scale);
            }
        }
        w[i] = (float)(scale * g);
    
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) 
        {
            for (k = l; k < n; k++) 
                scale += fabs((double)a[i][k]);
            if (scale) 
            {
                for (k = l; k < n; k++) 
                {
                    a[i][k] = (float)((double)a[i][k]/scale);
                    s += ((double)a[i][k] * (double)a[i][k]);
                }
                f = (double)a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (float)(f - g);
                for (k = l; k < n; k++) 
                    rv1[k] = (double)a[i][k] / h;
                if (i != m - 1) 
                {
                    for (j = l; j < m; j++) 
                    {
                        for (s = 0.0, k = l; k < n; k++) 
                            s += ((double)a[j][k] * (double)a[i][k]);
                        for (k = l; k < n; k++) 
                            a[j][k] += (float)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++) 
                    a[i][k] = (float)((double)a[i][k]*scale);
            }
        }
        anorm = MAX(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
    }
  
    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        if (i < n - 1) 
        {
            if (g) 
            {
                for (j = l; j < n; j++)
                    v[j][i] = (float)(((double)a[i][j] / (double)a[i][l]) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < n; k++) 
                        s += ((double)a[i][k] * (double)v[k][j]);
                    for (k = l; k < n; k++) 
                        v[k][j] += (float)(s * (double)v[k][i]);
                }
            }
            for (j = l; j < n; j++) 
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    

    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        l = i + 1;
        g = (double)w[i];
        if (i < n - 1) 
            for (j = l; j < n; j++) 
                a[i][j] = 0.0;
        if (g) 
        {
            g = 1.0 / g;
            if (i != n - 1) 
            {
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < m; k++) 
                        s += ((double)a[k][i] * (double)a[k][j]);
                    f = (s / (double)a[i][i]) * g;
                    for (k = i; k < m; k++) 
                        a[k][j] += (float)(f * (double)a[k][i]);
                }
            }
            for (j = i; j < m; j++) 
                a[j][i] = (float)((double)a[j][i]*g);
        }
        else 
        {
            for (j = i; j < m; j++) 
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) 
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++) 
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) 
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm) 
                {
                    flag = 0;
                    break;
                }
                if (fabs((double)w[nm]) + anorm == anorm) 
                    break;
            }
            if (flag) 
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) 
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) 
                    {
                        g = (double)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (float)h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) 
                        {
                            y = (double)a[j][nm];
                            z = (double)a[j][i];
                            a[j][nm] = (float)(y * c + z * s);
                            a[j][i] = (float)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)w[k];
            if (l == k) 
            {                  /* convergence */
                if (z < 0.0) 
                {              /* make singular value nonnegative */
                    w[k] = (float)(-z);
                    for (j = 0; j < n; j++) 
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) 
            {
                free((void*) rv1);
                fprintf(stderr, "No convergence after 30,000! iterations \n");
                return(0);
            }
    
            /* shift from bottom 2 x 2 minor */
            x = (double)w[l];
            nm = k - 1;
            y = (double)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
          
            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) 
            {
                i = j + 1;
                g = rv1[i];
                y = (double)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) 
                {
                    x = (double)v[jj][j];
                    z = (double)v[jj][i];
                    v[jj][j] = (float)(x * c + z * s);
                    v[jj][i] = (float)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = (float)z;
                if (z) 
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) 
                {
                    y = (double)a[jj][j];
                    z = (double)a[jj][i];
                    a[jj][j] = (float)(y * c + z * s);
                    a[jj][i] = (float)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (float)x;
        }
    }
    free((void*) rv1);
    return(1);
}
