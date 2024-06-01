#include <fstream>
#include <iostream>
using namespace std;

class morphology {
public:
    int numImgR;
    int numImgC;
    int imgMin;
    int imgMax;
    int numStructR;
    int numStructC;
    int structMin;
    int structMax;
    int rowOri;
    int colOri;

    int rFSize;
    int cFSize;
    int extR;
    int extC;
    int rSize;
    int cSize;

    int** ZFAry;
    int** morphAry;
    int** tempAry;
    int** structAry;

    morphology(ifstream& in, ifstream& struc) {
        in >> numImgR >> numImgC >> imgMin >> imgMax;
        struc >> numStructR >> numStructC >> structMin >> structMax >>
            rowOri >> colOri;

        rFSize = numStructR / 2;
        cFSize = numStructC / 2;
        extR = rFSize * 2;
        extC = cFSize * 2;
        rSize = numImgR + extR;
        cSize = numImgC + extC;

        ZFAry = new int* [rSize];
        morphAry = new int* [rSize];
        tempAry = new int* [rSize];
        structAry = new int* [numStructR];
        for (int i = 0; i < rSize; i++) {
            ZFAry[i] = new int[cSize];
            morphAry[i] = new int[cSize];
            tempAry[i] = new int[cSize];
        }
        for (int i = 0; i < numStructR; i++) {
            structAry[i] = new int[numStructC];
        }
    }

    void zero2DAry(int** ary, int r, int c) {
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                ary[i][j] = 0;
            }
        }
    }

    void loadImg(ifstream& in) {
        for (int i = rowOri; i < rSize - rFSize; i++) {
            for (int j = colOri; j < cSize - cFSize; j++) {
                in >> ZFAry[i][j];
            }
        }
    }

    void loadStruct(ifstream& s) {
        for (int i = 0; i < numStructR; i++) {
            for (int j = 0; j < numStructC; j++) {
                s >> structAry[i][j];
            }
        }
    }

    void computeDilation(int** inAry, int** outAry) {
        int i = rFSize;
        while (i < rSize) {
            int j = cFSize;
            while (j < cSize) {
                if (inAry[i][j] > 0) {
                    onePixelDilation(i, j, inAry, outAry);
                }
                j++;
            }
            i++;
        }
    }
    void computeErosion(int** inAry, int** outAry) {
        int i = rFSize;
        while (i < rSize) {
            int j = cFSize;
            while (j < cSize) {
                if (inAry[i][j] > 0) {
                    onePixelErosion(i, j, inAry, outAry);
                }
                j++;
            }
            i++;
        }
    }

    void onePixelDilation(int i, int j, int** inAry, int** outAry) {
        int iOffset = i - rowOri;
        int jOffset = j - colOri;
        int rIndex = 0;
        while (rIndex < numStructR) {
            int cIndex = 0;

            while (cIndex < numStructC) {
                if (structAry[rIndex][cIndex] > 0) {
                    outAry[iOffset + rIndex][jOffset + cIndex] = 1;
                }
                cIndex++;
            }
            rIndex++;
        }
    }

    void onePixelErosion(int i, int j, int** inAry, int** outAry) {
        int iOffset = i - rowOri;
        int jOffset = j - colOri;
        bool matchFlag = true;
        int rIndex = 0;
        while (matchFlag && rIndex < numStructR) {
            int cIndex = 0;

            while (matchFlag && cIndex < numStructC) {

                if (structAry[rIndex][cIndex] > 0 &&
                    inAry[iOffset + rIndex][jOffset + cIndex] <= 0) {
                    matchFlag = false;
                }
                cIndex++;
            }
            rIndex++;
        }
        if (matchFlag) {
            outAry[i][j] = 1;
        }
        else {
            outAry[i][j] = 0;
        }
    }

    void computeClosing(int** zAry, int** mAry, int** tAry) {
        computeDilation(zAry, tAry);
        computeErosion(tAry, mAry);
    }

    void computeOpening(int** zAry, int** mAry, int** tAry) {
        computeErosion(zAry, tAry);
        computeDilation(tAry, mAry);
    }

    void AryToFile(int** ary, ofstream& out) {
        out << numImgR << " " << numImgC << " " << imgMin << " " << imgMax << endl;
        for (int i = rFSize; i < rSize - rFSize; i++) {
            for (int j = cFSize; j < cSize - cFSize; j++) {
                out << ary[i][j] << " ";
            }
            out << endl;
        }
    }

    void prettyPrint(int** ary, ofstream& out) {
        out << numImgR << " " << numImgC << " " << imgMin << " " << imgMax << endl;
        for (int i = rFSize; i < rSize - rFSize; i++) {
            for (int j = cFSize; j < cSize - cFSize; j++) {
                if (ary[i][j] == 0) {
                    out << ". ";
                }
                else {
                    out << "1 ";
                }
            }
            out << endl;
        }
    }

    void prettyPrintStru(int** ary, ofstream& out) {
        out << numStructR << " " << numStructC << " " << structMin << " " << structMax << endl;
        for (int i = 0; i < numStructR; i++) {
            for (int j = 0; j < numStructC; j++) {
                if (ary[i][j] == 0) {
                    out << ". ";
                }
                else {
                    out << "1 ";
                }
            }
            out << endl;
        }
    }
};
int main(int argc, char* argv[]) {
    ifstream in(argv[1]);
    ifstream struc(argv[2]);
    ofstream dilateOutFile(argv[3]);
    ofstream erodeOutFile(argv[4]);
    ofstream openingOutFile(argv[5]);
    ofstream closingOutFile(argv[6]);
    ofstream prettyPrintFile(argv[7]);
    morphology m = morphology(in, struc);

    m.zero2DAry(m.ZFAry, m.rSize, m.cSize);
    m.loadImg(in);
    prettyPrintFile << "Printing Zero Framed Array" << endl;
    m.prettyPrint(m.ZFAry, prettyPrintFile);

    m.zero2DAry(m.structAry, m.numStructR, m.numStructC);
    m.loadStruct(struc);
    prettyPrintFile << "Printing structAry" << endl;
    m.prettyPrintStru(m.structAry, prettyPrintFile);

    m.zero2DAry(m.morphAry, m.rSize, m.cSize);
    m.computeDilation(m.ZFAry, m.morphAry);
    m.AryToFile(m.morphAry, dilateOutFile);
    prettyPrintFile << "Printing morphAry After Dilation" << endl;
    m.prettyPrint(m.morphAry, prettyPrintFile);

    m.zero2DAry(m.morphAry, m.rSize, m.cSize);
    m.computeErosion(m.ZFAry, m.morphAry);
    m.AryToFile(m.morphAry, erodeOutFile);
    prettyPrintFile << "Printing morphAry After Erosion" << endl;
    m.prettyPrint(m.morphAry, prettyPrintFile);

    m.zero2DAry(m.morphAry, m.rSize, m.cSize);
    m.computeOpening(m.ZFAry, m.morphAry, m.tempAry);
    m.AryToFile(m.morphAry, openingOutFile);
    prettyPrintFile << "Printing morphAry After Opening" << endl;
    m.prettyPrint(m.morphAry, prettyPrintFile);

    m.zero2DAry(m.morphAry, m.rSize, m.cSize);
    m.computeClosing(m.ZFAry, m.morphAry, m.tempAry);
    m.AryToFile(m.morphAry, closingOutFile);
    prettyPrintFile << "Printing morphAry After Closing" << endl;
    m.prettyPrint(m.morphAry, prettyPrintFile);


    prettyPrintFile.close();
    closingOutFile.close();
    erodeOutFile.close();
    openingOutFile.close();
    dilateOutFile.close();
    struc.close();
    in.close();
}