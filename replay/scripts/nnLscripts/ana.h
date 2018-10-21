#ifndef ana_h
#define ana_h 1

using namespace std;

#include "Setting.h"
#include "Tree.h"

/* +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+ */
class ana
{
  public:
    ana();
    ~ana();

  private:
    Setting *set;
    Tree *tr;

    int GetMaxEvent()     { return ENumMax; }

  public:
    void ReadTree(string name);
    void Roop();
    void Draw();

    bool Close();

    void MakeHist();

    void SetMaxEvent( int N )  { ENumMax = N; }

  private:
    int ENumMax;
    int ENum;

    bool GetEntry(int n);

    int LineID;
    int TextID;
    void DrawLine(double xmin, double ymin, double xmax, double ymax, int color=2, int style=1);
    void DrawText(double x, double y, string str, double size=0.05, int color=1);

  private:
// Lines, Textx
    TLine  *line[2000];
    TLatex *text[2000];

};

#endif
