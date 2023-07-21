#include <iostream>
#include <conio.h>
#include<stdlib.h>
#include<dos.h>
#include<math.h>
#include<iomanip>
unsigned int i,j,k,l,ch;

void banner(){
	cout<<endl<<"####################";
    cout<<endl<<"###              ###";
    cout<<endl<<"###    Welcome   ###";
    cout<<endl<<"###              ###";
    cout<<endl<<"####################";

    sleep(1);cout<<endl;
    cout<<endl<<"*** Math Problem ***";
    cout<<endl<<"*****  Solver  *****";
}

class stat{
    public:
    float x[100],y[100];
    unsigned int n,m;
    float sumx[5],sumy[5],sumxy[5];
    void getdata(){
      cout<<"Enter number of entries:";
      cin>>n;
      cout<<"Enter table of x :\n";
      for(i=0;i<n;i++)
        cin>>x[i];
      cout<<"Enter table of y :\n";
      for(i=0;i<n;i++)
        cin>>y[i];
    }
    stat(){
        for(i=0;i<5;i++){
            sumx[i]=0;
            sumy[i]=0;
            sumxy[i]=0;
        }
        m=5;
    }
    void sum();
    void line();
    float bxy();
    float byx();
    void regression();
    void corelation();
    void putdata(){
    	cout<<endl;
    	cout<<"x";
        for(i=0;i<n;i++)
        cout<<"\t"<<x[i];
    	cout<<endl<<"y";
        for(i=0;i<n;i++)
        cout<<"\t"<<y[i];
    }
    void putsum(){
    	cout<<endl;
        for(i=0;i<m;i++){
                cout<<endl<<"Sum of (x^"<<(i+1)<<")="<<sumx[i]<<" "
                <<endl<<"Sum of (y^"<<(i+1)<<")="<<sumy[i]<<" "
                <<endl<<"Sum of ((xy)^"<<(i+1)<<")="<<sumxy[i]<<endl;
        }
    }
};

class matrix{
    public:
    int mat[5][5];
    unsigned int m,n;
    void getsize(){
      cout<<"Enter size of matrix:"<<endl;
      cout<<"Number of rows:";
      cin>>m;
      cout<<"Number of columns:";
      cin>>n;
    }
    void getdata(){
      cout<<"Enter "<<(m*n)<<" elements :\n\n";
      for(i=0;i<m;i++)
          for(j=0;j<n;j++)
              cin>>mat[i][j];
    }
    void putdata(){
    	cout<<endl;
    	for(i=0;i<m;i++){
        	for(j=0;j<n;j++)
            	cout<<mat[i][j]<<" ";
        	cout<<endl;
        }
    }
}mtx1,mtx2,sum;

void stat::sum(){
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            sumx[i]=sumx[i]+pow(x[j],i+1);
            sumy[i]=sumy[i]+pow(y[j],i+1);
            sumxy[i]=sumxy[i]+pow(x[j]*y[j],(i+1));
        }
    }
}


void stat::line(){
    float a,b;
    m=2;
    sum();
    b=((sumy[0]*sumx[0]) - (sumxy[0]*n))/((sumx[0]*sumx[0] )- (sumx[1]*n));
    a=(sumy[0] - b*sumx[0])/n;
    setprecision(3);
    cout<<"y = "<<a<<" + "<<b<<"x";
}
//statestics

void num_integ(stat s){
    float x,r=0,h;
    h=s.x[1]-s.x[0];
    x=s.y[0]+s.y[s.n-1];
    for(i=1;i<s.n-1;i++)
        r=r+s.y[i];
    cout<<((x+2*r)*h)/2;
}
    
float stat::bxy(){
    setprecision(3);
    return ((sumxy[0] / n) -  (sumx[0]*sumy[0])/(n*n)) / ((sumy[1]/ n) - (sumy[0]*sumy[0])/ (n*n));
}

float stat::byx(){
    setprecision(3);
    return ((sumxy[0] / n) - (sumx[0]*sumy[0])/(n*n)) / ((sumx[1]/ n) - (sumx[0]*sumx[0])/ (n*n));
}

void stat::regression(){
    cout<<endl;
    cout<<"Line of regression of y on x :-"<<endl;
    setprecision(3);
    cout<<"y = "<<byx()<<" x + "<<((-(sumx[0]/n)*byx()) + (sumy[0]/n));
    
    cout<<endl<<"Line of regression of x on y :-"<<endl;  
    setprecision(3);
    cout<<"x = "<<bxy()<<" y - "<<((-(sumy[0]/n)*bxy()) + (sumx[0]/n));
}

void stat::corelation(){        
    setprecision(3);
    cout<<endl<<"Coefficient of corelation = "<<sqrt(bxy()*byx());
}

//matrix

void echelon(){    
    unsigned int m1,m2;
    for(k=0;k<mtx1.n && k<mtx1.m;k++)
        for(j=mtx1.m;j>k;j--){
            m1=mtx1.mat[k][k];      
            m2=mtx1.mat[j][k];
            if(mtx1.mat[k][k]==mtx1.mat[j][k])
                for(i=0;i<mtx1.n;i++)
                    mtx1.mat[j][i]=mtx1.mat[j][i]-mtx1.mat[k][i];
            else    
                for(i=0;i<mtx1.n;i++)
                    mtx1.mat[j][i]=mtx1.mat[j][i]*m1-m2*mtx1.mat[k][i];
    }
    unsigned int count1=0,count2=0;
    for(i=0;i<mtx1.m;i++){
        count1=0;
        for(j=0;j<mtx1.n;j++)
            if(mtx1.mat[i][j]==0)
                count1++;
        if(count1==mtx1.n)
            count2++;
    }
    cout<<(mtx1.m-count2);
}

float determinant(matrix mtx){
    if(mtx.m==mtx.n){
        float det=0;
        sum.m=mtx.m-1;
        if(mtx.m==2)
            return mtx.mat[0][0]*mtx.mat[1][1]-mtx.mat[1][0]*mtx.mat[0][1];
        else
        for(i=0;i<mtx.m;i++){
            int temp1=0;
            for(j=1;j<mtx.m;j++){
                int temp2=0;
                for(k=0;k<mtx.m;k++){
                    if(i==k)
                        continue;
                    sum.mat[temp1][temp2]=mtx.mat[j][k];
                    temp2++;
                }
                temp1++;
            }
            det=det+(pow(-1,i)*mtx.mat[0][i]*determinant(sum));
        }
        return det;
    }
    else{
        cout<<"Size of matrix is wrong. (Note :- Matrix has to be square matrix)";
        return 0;
    }
}

matrix minor(matrix mtx){
    if(mtx.m==mtx.n){
        matrix mnr;
        float min;
        mnr=mtx;
        sum.m=mtx.m-1;
        if(mtx.m==2){
            min=mtx.mat[0][0];
            mtx.mat[0][0]=mtx.mat[1][1];
            mtx.mat[1][1]=min;
            mtx.mat[1][0]=-mtx.mat[1][0];
            mtx.mat[0][1]=-mtx.mat[0][1];
            mtx.putdata();
        }
        else
        for(l=0;l<mtx.m;l++)
            for(i=0;i<mtx.m;i++){
                int temp1=0;
                for(j=0;j<mtx.m;j++){
                    int temp2=0;
                    for(k=0;k<mtx.m;k++){
                        if(i==k || l==j)
                            continue;
                        sum.mat[temp1][temp2]=mtx.mat[j][k];
                        temp2++;
                    }
                    if(l==j)
                        continue;
                    temp1++;
                }
                mnr.mat[l][i]=determinant(sum);
            }
        return mnr;
    }
    else{
        cout<<"Size of matrix is wrong. (Note :- Matrix has to be square matrix)";
        return mtx;
    }
}

float root1(float a,float b,float c){
    return ((-b+sqrt(b*b-4*a*c))/(2*a));
}
float root2(float a,float b,float c){
    return ((-b-sqrt(b*b-4*a*c))/(2*a));
}

//calculation

//Addition of matrix

void add(){
	sum.m=mtx1.m;
    sum.n=mtx1.n;
    for (i=0;i<mtx1.m;i++ )
    	for (j=0;j<mtx1.n;j++)
        	sum.mat[i][j]=mtx1.mat[i][j]+mtx2.mat[i][j];
    sum.putdata();
}

//substraction of matrix

void sub(){
	sum.m=mtx1.m;
    sum.n=mtx1.n;
    for (i=0;i<mtx1.m;i++ )
    	for (j=0;j<mtx1.n;j++)
        	sum.mat[i][j]=mtx1.mat[i][j]-mtx2.mat[i][j];
    sum.putdata();
}

//Multiplication of matrix
void multi(){
  if(mtx1.n==mtx2.m){
        sum.m=mtx1.m;
        sum.n=mtx2.n;
        for(i=0;i<mtx1.m;i++)
            for(j=0;j<mtx2.n;j++){
                sum.mat[i][j]=0;
                for(k=0;k<mtx1.n;k++)
                    sum.mat[i][j]+=mtx1.mat[i][k]*mtx2.mat[k][j];
            }
        sum.putdata();
    }
    else
        cout<<"Size of matrix is wrong.  (Note :- Number of columns of first matrix and Number of rows of second matrix must be same.)";
}

//eigen value
void eigen(){
    matrix m;
    int n1;float temp[4],lamda[3];
        if(mtx1.m==mtx1.n){
            if(mtx1.m==2){
                temp[0]=1;
                temp[1]=-mtx1.mat[0][0]-mtx1.mat[1][1];
                temp[2]=determinant(mtx1);
                }
            if(mtx1.m==3){
                int s[4]={1,0};
                for(i=0;i<mtx1.m;i++)
                    s[1]=s[1]+mtx1.mat[i][i];
                s[1]=-s[1];
                m=minor(mtx1);
                for(i=0;i<m.m;i++)
                    s[2]+=m.mat[i][i];
                s[3]=-determinant(mtx1);
                if(s[3]<0)
                    n1=-s[3];
                else
                    n1=s[3];
                for(int i=-n1;i<n1;i++){
                    temp[0]=1;
                    if(i==0)
                      continue;
                    for(j=1;j<4;j++)
                        temp[j]=temp[j-1]*i+s[j];
                    if(temp[3]==0){
                        lamda[0]=i;
                        cout<<lamda[1]<<",";
                        break;
                    }
                }
            }
            lamda[1]=root1(temp[0],temp[1],temp[2]);
            lamda[2]=root2(temp[0],temp[1],temp[2]);
            for(i=1;i<3;i++){
                cout<<lamda[i];
                if(i<2)
                    cout<<",";
        }
    }
}

void calculator(){
    float a,b;
    char t;
    char ch;
    do{
    jump2:
        cout<<"-------------------------------";
        cout<<"\nCalculate :\n";
        cin>>a;
        int a1=a,b1;
        cin>>ch;
        switch(ch){
            case '+':
    	        cin>>b;
        	    a=a+b;
            	break;

            case '-':
	            cin>>b;
    	        a=a-b;
        	    break;

            case '*':
	            cin>>b;
            	a=a*b;
        	    break;

            case '/':
    	        cin>>b;
	            a=a/b;
            	break;

            case '%':
            	cin>>b1;
                a=a1%b1;
                break;

            case '^':
                cin>>b;
                a=pow(a,b);
                break;

            case '!':
	            if(a>=0){
    	            if(a==0){
        	            a=1;
            	        break;
                	}
                	for(i=a-1;i>1;i--)
                	    a=a*i;
            	}
            	else
                	cout<<"Factorial of a negative number can not be calculated.";
	            break;

            default:
            cout<<"Your choise is wrong.Try again:";
            system("cls");
            goto jump2;
        }
        cout<<'='<<a<<endl<<"Want to continue calculation?(y/n)";
        cin>>t;
    }while(t=='y' || t=='Y');
}

void matrix_opt(){
    int ch;                                                                                                   
    jump:
    cout<<"Matrix :-"<<endl<<endl;
    cout<<"1.Addition\n2.Substraction\n3.Multiplication\n4.Determinant\n5.Minor\n6.Rank of matrix\n7.Eigen Value\n8.Back"<<endl<<endl;
    cout<<"Enter your choice:";
    cin>>ch;
    system("cls");
    if(ch<8 && ch>0){
	    mtx1.getsize();
    	mtx1.getdata();
    }
    switch(ch)
    {
        case 1:
            mtx2=mtx1;
            cout<<"Another matrix:-"<<endl;
            mtx2.getdata();           
            system("cls");
            cout<<"Addition of matrix "<<endl;
            mtx1.putdata();
            cout<<endl<<"and"<<endl;
            mtx2.putdata();
            cout<<endl<<"is"<<endl;
            add();
            getch();
            break;

        case 2:
            mtx2=mtx1;
            cout<<"Another matrix:-"<<endl;
            mtx2.getdata();
            system("cls");
            cout<<"Substraction of matrix "<<endl;
            mtx1.putdata();
            cout<<endl<<"and"<<endl;
            mtx2.putdata();
            cout<<endl<<"is"<<endl;
            sub();
            getch();
            break;

        case 3:
            cout<<"Another matrix:-"<<endl;
            mtx2.getsize();
            mtx2.getdata();
            system("cls");
            cout<<"Multiplication of matrix "<<endl;
            mtx1.putdata();
            cout<<endl<<"and"<<endl;
            mtx2.putdata();
            cout<<endl<<"is"<<endl;
            multi();
            getch();
            break;

        case 4:         
            system("cls");
            cout<<"Determinant of matrix";
            mtx1.putdata();
            cout<<endl<<" is "<<determinant(mtx1);
            getch();
            break;

        case 5:
            mtx1=minor(mtx1);
            system("cls");
            cout<<"Minors of matrix"<<endl;
            mtx1.putdata();
            cout<<endl<<" are ";
            mtx1.putdata();
            getch();
            break;

        case 6:              
            system("cls");
        	cout<<"Rank of matrix is";
            mtx1.putdata();
            cout<<endl<<" is ";
            echelon();
            getch();
            break;

        case 7:                
            system("cls");
            cout<<"Eigen values of matrix"; 
            mtx1.putdata();
            cout<<endl<<" are ";
            eigen();                  
            getch();
            break;

        case 8:
        	break;

        default:
            cout<<"You choice is wrong.Try again"<<endl;
            goto jump;
    }
}

void stat_opt(){
    stat s;
    jump:
    system("cls");
    cout<<"Statistics:-"<<endl<<endl;
    cout<<"1.sum of (x , y , xy)\n2.Fit straight line\n3.Line of regression\n4.Coefficient of corelation\n5.Numerical integration\n6.Back";
    cout<<endl<<endl<<"Select above option:";
    cin>>ch;
    system("cls");
    if(ch<6 && ch>0){
	    s.getdata();
    	s.sum();
    }
    switch(ch)
    {
        case 1:
            cout<<"Enter power upto which you want sumation:-";
            cin>>s.m;
            system("cls");
            s.putdata();
            s.putsum();
            getch();
            break;

        case 2:                
            system("cls");
            s.putdata();
            cout<<endl<<"Equation of line is"<<endl;
            s.line();
            getch();
            break;

        case 3:                
            system("cls");
            s.putdata();
            s.regression();
            getch();
            break;

        case 4:                
            system("cls");
            s.putdata();
            s.corelation();
            getch();
            break;

        case 5:                
            system("cls");
            s.putdata();
            cout<<endl<<"Numerical integration is ";
            num_integ(s);
            getch();
            break;

        case 6:
        	break;

        default:
            cout<<"You choice is wrong.Try again"<<endl;
            goto jump;
    }
}

int main()
{
	banner();
    sleep(3);
    do{
    	system("cls");
        cout<<"\n1.Simple Calculator\n2.Matrix\n3.Statistics\n4.exit\n";
        cout<<"Select above option:";
        cin>>ch;
        sleep(1);
        system("cls");
        switch(ch){
            case 1:
                calculator();
                break;

            case 2:
                matrix_opt();
                break;

            case 3:
                stat_opt();
                break;

            case 4:
                cout<<"Exiting";
                for(i=0;i<5;i++)
                {
                	sleep(1);
                    cout<<".";
                }
                break;

            default:
                cout<<"You choice is wrong.Try again";
                sleep(3);
        }
    }while(ch!=4);
    return 0;
}
