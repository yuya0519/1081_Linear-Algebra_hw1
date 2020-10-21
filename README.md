# 1081_Linear-Algebra_hw1 
## 主要目標:  
分析 Protein-DNA 鍵結位置(Protein-DNA binding sites)。

## 程式運行方式:  

**該程式執行軟體為Visual 2019。**   
請將程式碼複製至上述軟體中，再編譯執行。  
創建3個檔案，分別為 1.DNA存放DNA結構、1.PROTEIN存放蛋白質結構和1.MATEIX存放要用來轉換的平移矩陣並將3個檔案放入和程式同一個資料夾。

## 簡要使用說明:  
- Input: 輸入皆由讀檔，檔案格式為: 
- 1.DNA  
![avatar](https://upload.cc/i1/2020/10/21/zfPIZu.jpg) 

- 1.PROTEIN  
![avatar](https://upload.cc/i1/2020/10/21/vNUTmF.jpg) 

- 1.MATRIX  
![avatar](https://upload.cc/i1/2020/10/21/sBqXZb.jpg) 


- Output: 第一行為鍵結個數，第二行為鍵結編號且每個位置編號皆空一格。 
> 輸出範例  
![avatar](https://upload.cc/i1/2020/10/21/09Ie12.jpg)  

## 程式結構說明:  

### 基本原理:  
利用 Protein 結構、DNA 結構以及平移旋轉矩陣將 protein 轉成適當位置，並計算 Protein 與 DNA 之間的距離，若 Protein 與 DNA 之間的距離小於 4.5Å，則此胺基酸(aminoacid)位置為 Protein-DNA鍵結位置。 

### 引用函式庫說明:
```cpp
#include<iostream> //負責輸出/輸入
#include<string> //負責處理string相關函式
#include<fstream> //負責讀檔
#include<vector> //引用vector
#include<cmath> //為引用pow()以計算次方、sqrt()用以開根號
#include<iomanip> //用來排版
#include<cstdlib> //將相關名稱新增至 std 命名空間
```

### Global變數說明:
```cpp
stream file2("1.DNA");  //讀檔 
ifstream file1("1.PROTEIN");    //讀檔 
ifstream file3("1.MATRIX"); //讀檔 
string world;   //用以存放讀入的字串 
char w[20]; //用來暫存分割字串 

vector<double>buf;  //暫存一組的分割字串
vector<vector<double> > dnaPosition;    //存放所有DNA資料 
vector<vector<double> > proPosition;    //存放所有protein資料
vector<vector<double> > matrix; //存放平移矩陣
vector<double>proPosition_x;    //存放protein x軸位置
vector<double>proPosition_y;    //存放protein y軸位置
vector<double>proPosition_z;    //存放protein z軸位置
vector<int>proPosition_aa;  //存放已檢查的protein編號
vector<int>  binding_sites; //存放有鍵結的protein編號
```
### 詳細程式碼說明

```cpp
while (getline(file1, world)) {
		for (int i = 22; i < 54; i++) {
			if (world[i] != '\0') {
				int m = 0;
				while (world[i] != '\0') {
					w[m] = world[i];
					m++;
					if (i == 26 || i == 38 || i == 46 || i == 54)
						break;
					i++;
				}
				buf.push_back(atof(w));
				for (int i = 0; i < 20; i++)
					w[i] = ' ';
			}
		}
		proPosition.push_back(buf);
		buf.clear();
	}
```
> 將讀入的porein資料做處理，並存入proPosition中。

```cpp
while (getline(file2, world)) {
		for (int i = 22; i < 54; i++) {
			if (world[i] != '\0') {			
				int m = 0;
				while (world[i] != '\0' ) {
					w[m] = world[i];
					m++;
					if (i == 26 || i == 38 || i == 46 || i == 54 )
						break;
					i++;
				}
				buf.push_back(atof(w));
				for (int i = 0; i < 20; i++)
					w[i] = ' ';
			}					
		}
		dnaPosition.push_back(buf);
		buf.clear();
	}
```
> 將讀入的DNA資料做處理，並存入dnaPosition中。
```cpp
    int e = 0;
	while (getline(file3, world)) {	
			if(e<6&&2<e){
				for (int i = 2; i < world.length(); i++) {
					if (world[i] != '\0'&& world[i]!='\t' && world[i] != ' ') {
						int m = 0;
						while (world[i] != '\0' && world[i] != '\t' && world[i] != '\n' && world[i] != ' ') {
							w[m] = world[i];
							m++;
							i++;
						}
						buf.push_back(atof(w));
						for (int i = 0; i < 20; i++)
							w[i] = ' ';
					}
				}
				matrix.push_back(buf);
				buf.clear();
		}
			e++;			
	}
```
> 將讀入的matrix資料做處理，並存入matrix中。
```cpp
	for (int i = 0; i < proPosition.size(); i++) {
		proPosition_aa.push_back(proPosition[i][0]);
		proPosition_x.push_back(matrix[0][0] + (matrix[0][1] * proPosition[i][1]) + (matrix[0][2] * proPosition[i][2]) + (matrix[0][3] * proPosition[i][3]));
		proPosition_y.push_back( matrix[1][0] + (matrix[1][1] * proPosition[i][1]) + (matrix[1][2] * proPosition[i][2]) + (matrix[1][3] * proPosition[i][3]));
		proPosition_z.push_back( matrix[2][0] + (matrix[2][1] * proPosition[i][1]) + (matrix[2][2] * proPosition[i][2]) + (matrix[2][3] * proPosition[i][3]));
		for (int j = 0; j < dnaPosition.size(); j++) {
					if ((sqrt(pow((proPosition_x[i] - dnaPosition[j][1]), 2) + pow((proPosition_y[i] - dnaPosition[j][2]), 2) + pow((proPosition_z[i] - dnaPosition[j][3]), 2))) < 4.5) {
						if (binding_sites.size() == 0) {
							binding_sites.push_back(proPosition_aa[i]);
						}
						else if (binding_sites[binding_sites.size() - 1] != proPosition_aa[i]) {
							binding_sites.push_back(proPosition_aa[i]);
						}
						break;
					}
				}
	}
```
> 利用計算距離的公式計算，
> 並判斷該距離是否符合鍵結標準
> 存入binding_sites中
```cpp
	cout << binding_sites.size() << endl;
	for (int i = 0; i < binding_sites.size(); i++) {
		if (i != binding_sites.size() - 1)
			cout << binding_sites[i] << " ";
		else
			cout << binding_sites[i];
	}
	cout << endl;
	file1.close();
	file2.close();
	file3.close();
```
> 依照格式輸出