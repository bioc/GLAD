
struct paire_double
{
  double value;
  int index_add;
};

struct split_region
{
  vector<double> LogRatio;
  vector<int> index;
};



struct agg
{
  double Mean;
  double Var;
  double VarLike;
  int Card; 
  int LabelRegion;

  agg() {
    Mean = Var = VarLike = 0;
    Card = 0;
  }

  agg(double Mean, double Var, double VarLike, int Card, int LabelRegion) :
    Mean(Mean), Var(Var), VarLike(VarLike), Card(Card), LabelRegion(LabelRegion)
  {
  }

#ifndef USE_COMPARE
  bool operator<(const agg &b) const {
    return LabelRegion < b.LabelRegion;
  }
#endif
};
