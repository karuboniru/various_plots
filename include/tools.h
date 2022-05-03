#pragma once
#include <TAxis.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPaletteAxis.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <string>

namespace style
{
    static constexpr Double_t fgkTextSize = 0.05;
    static constexpr Double_t fgkTitleSize = 0.05;
    static constexpr Double_t fgkMarkerSize = 1;
    static constexpr Double_t fgkLineWidth = 2;
    static constexpr Int_t fgkTextFont = 42;
    static constexpr Double_t fgkLabelOffset = 0.01;
    static constexpr Double_t fgkXTitleOffset = 1.25; // 1.1;//1.25;
    static constexpr Double_t fgkYTitleOffset = 1.1;  // 1.2;
    static constexpr Double_t fgkTickLength = 0.02;
} // namespace name

class global_style
{
public:
    global_style()
    {
        using namespace style;
        std::cout << "Setting Style" << std::endl;
        gStyle->SetFrameBorderMode(0);
        gStyle->SetFrameFillColor(0);
        gStyle->SetCanvasBorderMode(0);
        gStyle->SetPadBorderMode(0);
        gStyle->SetPadColor(10);
        gStyle->SetCanvasColor(10);
        gStyle->SetTitleFillColor(10);
        gStyle->SetTitleBorderSize(-1);
        gStyle->SetStatColor(10);
        gStyle->SetStatBorderSize(-1);
        // gStyle->SetLegendBorderSize(-1);
        //
        gStyle->SetDrawBorder(0);
        gStyle->SetTextFont(fgkTextFont);
        gStyle->SetStatFont(fgkTextFont);
        gStyle->SetStatFontSize(fgkTextSize);
        gStyle->SetStatX(0.97);
        gStyle->SetStatY(0.98);
        gStyle->SetStatH(0.03);
        gStyle->SetStatW(0.3);
        gStyle->SetTickLength(fgkTickLength, "xy");
        gStyle->SetEndErrorSize(3);
        gStyle->SetLabelSize(fgkTextSize, "xyz");
        gStyle->SetLabelFont(fgkTextFont, "xyz");
        gStyle->SetLabelOffset(fgkLabelOffset, "xyz");
        gStyle->SetTitleFont(fgkTextFont, "xyz");
        gStyle->SetTitleFont(fgkTextFont, "");
        gStyle->SetTitleFontSize(fgkTitleSize);
        gStyle->SetTitleOffset(fgkXTitleOffset, "x");
        gStyle->SetTitleOffset(fgkYTitleOffset, "y");
        gStyle->SetTitleOffset(1.0, "z");
        gStyle->SetTitleSize(fgkTitleSize, "xyz");
        gStyle->SetTitleSize(fgkTitleSize, "");
        gStyle->SetMarkerSize(fgkMarkerSize);
        gStyle->SetPalette(1, 0);
        TGaxis::SetMaxDigits(3);
        gStyle->SetTitleBorderSize(-1);
        gROOT->ForceStyle();
    }
};

inline global_style global_style_instance;

template <typename T>
void PadSetup(T &currentPad, const Double_t currentLeft = 0.12, const Double_t currentTop = 0.09, const Double_t currentRight = 0.13, const Double_t currentBottom = 0.14)
{
    currentPad->SetTicks(1, 1);
    currentPad->SetLeftMargin(currentLeft);
    currentPad->SetTopMargin(currentTop);
    currentPad->SetRightMargin(currentRight);
    currentPad->SetBottomMargin(currentBottom);

    currentPad->SetFillColor(0); // this is the desired one!!!
}

template <typename T>
void AxisStyle(T *ax, Bool_t kcen)
{
    using namespace style;
    static_assert(std::is_base_of_v<TAxis, T> || std::is_base_of_v<TGaxis, T>, "AxisStyle: T must be a TAxis or TGaxis");
    ax->SetTickLength(fgkTickLength);

    ax->SetLabelFont(fgkTextFont);
    ax->SetLabelSize(fgkTextSize);
    ax->SetLabelOffset(fgkLabelOffset);

    ax->SetTitleFont(fgkTextFont);
    ax->SetTitleSize(fgkTitleSize);

    kcen = 1;
    ax->CenterTitle(kcen);

    ax->SetNdivisions(505);
    if (std::is_base_of_v<TGaxis, T>)
    {
        ax->SetTitleOffset(fgkXTitleOffset);
    }
}

template <typename T,
          std::enable_if_t<std::is_base_of_v<TH1, T> || std::is_base_of_v<TH2, T> || std::is_base_of_v<THStack, T>, int> = 0>
void ResetStyle(T *obj, TVirtualPad *cpad, Bool_t kcen)
{
    using namespace style;
    if (!obj)
    {
        printf("style::ResetStyle obj null!\n");
        exit(1);
    }

    AxisStyle(obj->GetXaxis(), kcen);
    AxisStyle(obj->GetYaxis(), kcen);

    obj->GetXaxis()->SetTitleOffset(fgkXTitleOffset);
    obj->GetYaxis()->SetTitleOffset(fgkYTitleOffset);
    if constexpr (!std::is_base_of_v<THStack, T>)
    {
        obj->SetMarkerSize(fgkMarkerSize);
        if (cpad)
        {
            TPaletteAxis *palette = (TPaletteAxis *)obj->GetListOfFunctions()->FindObject("palette");
            if (!palette)
            {
                printf("ResetStyle no palette!!\n");
                obj->GetListOfFunctions()->Print();
            }
            else
            {
                palette->SetX1NDC(1 - cpad->GetRightMargin() + 0.005);
                palette->SetX2NDC(1 - cpad->GetRightMargin() / 3 * 2);
                palette->SetY1NDC(cpad->GetBottomMargin());
                palette->SetY2NDC(1 - cpad->GetTopMargin());
                palette->SetLabelFont(fgkTextFont);
                palette->SetLabelSize(fgkTextSize);
                palette->SetLabelOffset(fgkLabelOffset);
            }
        }
    }
}

decltype(auto) inline getCanvas(const char *name = "")
{
    const double factor = 2;
    auto c = std::make_unique<TCanvas>(name, name, 800 * factor, 600 * factor);
    c->cd();
    return c;
}

template <typename T>
void plot(const T &obj, const std::string &path_prefix, const char *opt = "", TLegend *leg = nullptr)
{
    using std::string_literals::operator""s;
    // std::unique_ptr<TCanvas> c{new TCanvas{}};
    auto c = getCanvas();
    PadSetup(c);
    if (opt == "colz"s)
    {
        c->SetRightMargin(0.15);
    }
    if constexpr (std::is_base_of_v<THStack, std::remove_cvref_t<decltype(*obj)>>)
    {
        // THStack need to be drawn first, or Get?axis() will return nullptr
        obj->Draw(opt);
        if (!obj->GetXaxis())
        {
            return; // do nothing for empty THStack, instead of crashing
        }
    }
    ResetStyle(&*obj, c.get()->GetPad(0), 1);
    obj->Draw(opt);
    if (leg)
    {
        leg->Draw();
    }
    c->Update();
    c->WaitPrimitive();
    auto path_png = path_prefix + obj->GetName() + ".png"s;
    auto path_pdf = path_prefix + obj->GetName() + ".pdf"s;
    c->SaveAs(path_png.c_str());
    c->SaveAs(path_pdf.c_str());
}

template <typename T>
void plot(const T &obj, const std::string &path_prefix, double xsec, const char *opt = "")
{
    using std::string_literals::operator""s;
    auto c = getCanvas();
    // std::unique_ptr<TCanvas> c{new TCanvas{}};
    PadSetup(c);
    if (opt == "colz"s)
    {
        c->SetRightMargin(0.15);
    }
    if constexpr (std::is_base_of_v<THStack, std::remove_cvref_t<decltype(*obj)>>)
    {
        // THStack need to be drawn first, unless Get?axis() will return nullptr
        obj->Draw(opt);
    }
    ResetStyle(&*obj, c.get()->GetPad(0), 1);
    obj->Draw(opt);
    c->Update();
    // std::unique_ptr<TLatex> la{new TLatex{0.15, 0.85, ("#sigma = "s + std::to_string(xsec) + "#times 10^{-38} cm^{2}").c_str()}};
    // la->DrawLatex(0.15, 0.85, ("#sigma = "s + std::to_string(xsec) + "#times 10^{-38} cm^{2}").c_str());
    // la->Draw();
    // la->Paint();
    std::unique_ptr<TLegend> leg{new TLegend};
    leg->AddEntry((TObject *)0, ("#sigma = "s + std::to_string(xsec) + "#times 10^{-38} cm^{2}").c_str(), "p");
    leg->SetBorderSize(-1);
    leg->Draw();
    c->Update();
    c->WaitPrimitive();
    auto path_png = path_prefix + obj->GetName() + ".png"s;
    auto path_pdf = path_prefix + obj->GetName() + ".pdf"s;
    c->SaveAs(path_png.c_str());
    c->SaveAs(path_pdf.c_str());
}

template <typename T, typename Y>
void plot_with_normalized(const T &obj, const Y &hist1d, const std::string &path_prefix, const char *opt = "", const bool h_on_x = true)
{
    static_assert(std::is_base_of_v<TH2, std::remove_cvref_t<decltype(*obj)>>, "T must be derived from TH2");
    static_assert(std::is_base_of_v<TH1, std::remove_cvref_t<decltype(*hist1d)>>, "Y must be derived from TH1");
    using std::string_literals::operator""s;
    // std::unique_ptr<TCanvas> c{new TCanvas{}};
    auto c = getCanvas();
    PadSetup(c);
    double max = h_on_x ? obj->GetYaxis()->GetXmax() : obj->GetXaxis()->GetXmax();
    double min = h_on_x ? obj->GetYaxis()->GetXmin() : obj->GetXaxis()->GetXmin();
    hist1d->Scale(1. / hist1d->GetMaximum() * (max - min) * 0.9);
    std::size_t n_binsx = h_on_x ? obj->GetNbinsX() : obj->GetNbinsY();
    std::unique_ptr<double[]> x_center{new double[n_binsx]}, value{new double[n_binsx]};
    for (std::size_t i = 0; i < n_binsx; ++i)
    {
        x_center[i] = h_on_x ? obj->GetXaxis()->GetBinCenter(i + 1) : obj->GetYaxis()->GetBinCenter(i + 1);
        value[i] = hist1d->GetBinContent(i + 1) + min;
    }
    std::unique_ptr<TGraph> graph{h_on_x ? new TGraph(n_binsx, x_center.get(), value.get()) : new TGraph(n_binsx, value.get(), x_center.get())};
    graph->SetLineWidth(2);
    obj->Draw(opt);
    graph->Draw("same");
    // hist1d->Draw(h_on_x ? "same" : "same hbar");
    // c->Update();
    // c->WaitPrimitive();
    auto path_png = path_prefix + obj->GetName() + "_new.png"s;
    auto path_pdf = path_prefix + obj->GetName() + "_new.pdf"s;
    c->SaveAs(path_png.c_str());
    c->SaveAs(path_pdf.c_str());
}

template <typename T>
auto normalize_slice(const T &hist, bool on_axis_x = true)
{
    typedef std::remove_cvref_t<decltype(*hist)> hist_type;
    static_assert(std::is_base_of_v<TH2, hist_type>, "T must be derived from TH2");
    using std::string_literals::operator""s;
    auto hist_name_new = hist->GetName() + (on_axis_x ? "_norm_x"s : "_norm_y"s);
    auto hist_norm = reinterpret_cast<hist_type *>(hist->Clone(hist_name_new.c_str()));
    auto axis = on_axis_x ? hist->GetXaxis() : hist->GetYaxis();
    auto o_axis = on_axis_x ? hist->GetYaxis() : hist->GetXaxis();
    for (int i = 1; i < axis->GetNbins() + 1; i++)
    {
        double max = 0;
        for (int j = 1; j < o_axis->GetNbins() + 1; j++)
        {
            auto [x, y] = on_axis_x ? std::make_pair(i, j) : std::make_pair(j, i);
            max = std::max(max, hist->GetBinContent(x, y));
        }
        for (int j = 1; j < o_axis->GetNbins() + 1; j++)
        {
            auto [x, y] = on_axis_x ? std::make_pair(i, j) : std::make_pair(j, i);
            if (max != 0)
                hist_norm->SetBinContent(x, y, hist->GetBinContent(x, y) / max);
        }
    }
    return std::unique_ptr<hist_type>(hist_norm);
}

template <typename T>
void normalize_plot(const T &obj, const std::string &path_prefix, const char *opt = "", const bool h_on_x = true)
{
    std::unique_ptr<TH1> h1{h_on_x ? obj->ProjectionX() : obj->ProjectionY()};
    auto normalized_hist = normalize_slice(obj, h_on_x);
    plot_with_normalized(normalized_hist, h1, path_prefix, opt, h_on_x);
}