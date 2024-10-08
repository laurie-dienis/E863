#ifndef ActSRIM_h
#define ActSRIM_h

#include "TF1.h"
#include "TSpline.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

    class SRIM
    {
    public:
        using PtrSpline = std::unique_ptr<TSpline3>;
        using PtrFunc = std::unique_ptr<TF1>;

    private:
        std::vector<std::string> fKeys; //!< Store known tables
        // Energy->Range
        std::map<std::string, std::unique_ptr<TSpline3>> fSplinesDirect;
        std::map<std::string, std::unique_ptr<TF1>> fInterpolationsDirect;
        // Range->Energy
        std::map<std::string, std::unique_ptr<TSpline3>> fSplinesInverse;
        std::map<std::string, std::unique_ptr<TF1>> fInterpolationsInverse;
        // Energy->Stopping powers (nuclear + electronic)
        std::map<std::string, std::unique_ptr<TSpline3>> fSplinesStoppings;
        std::map<std::string, std::unique_ptr<TF1>> fStoppings;
        // Range->Longitudinal straggling
        std::map<std::string, std::unique_ptr<TSpline3>> fSplinesLongStrag;
        std::map<std::string, std::unique_ptr<TF1>> fLongStrag;
        // Range->Lateral straggling
        std::map<std::string, std::unique_ptr<TSpline3>> fSplinesLatStrag;
        std::map<std::string, std::unique_ptr<TF1>> fLatStrag;


    public:
        SRIM() = default;
        SRIM(const std::string& material, const std::string& file);

        void ReadTable(const std::string& key, const std::string& file);

        [[deprecated("Favour use of new SRIM::ReadTable(...) which does not require manual edition of SRIM file")]] void
        ReadInterpolations(std::string key, std::string fileName);

        double EvalDirect(const std::string& key, double energy) { return fInterpolationsDirect[key]->Eval(energy); }
        double EvalInverse(const std::string& key, double range) { return fInterpolationsInverse[key]->Eval(range); }

        // Evaluate the other columns of the SRIM table
        double EvalStoppingPower(const std::string& key, double energy) { return fStoppings[key]->Eval(energy); }
        double EvalLongStraggling(const std::string& key, double range) { return fLongStrag[key]->Eval(range); }
        double EvalLatStraggling(const std::string& key, double range) { return fLatStrag[key]->Eval(range); }

        void Draw(const std::string& what, const std::vector<std::string>& keys = {});

        double Slow(const std::string& material, double Tini, double thickness, double angleInRad = 0);

        double EvalInitialEnergy(const std::string& material, double Tafter, double thickness, double angleInRad = 0);

        bool CheckKeyIsStored(const std::string& key);

    private:
        bool IsBreakLine(const std::string& line);
        double ConvertToDouble(std::string& str, const std::string& unit);
        PtrSpline GetSpline(std::vector<double>& x, std::vector<double>& y, const std::string& name);
    };
#endif
