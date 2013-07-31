#include <phypp.hpp>
#include <plot.hpp>

int main(int argc, char* argv[]) {
    {
        plot::eps p("out/test_color.eps", 200, 200);

        vec1rgb c = {rgb(1,0,0), rgb(0,1,0), rgb(0,0,1), rgb(0,0,0)};

        p.line_width(2);
        p.font(12);
        p.line_join(plot::join::round);

        p.view_grid(2, 2);
        for (int_t i = 0; i < 4; ++i) {
            p.color(c[i]);
            p.start_path().move_to(10,10).line(20,0).line(0,30).line(-10,10).line(-10,-10).close_path();
            p.stroke();
            p.move_to(50,50).show("Hello\nyou", 0.5);
            p.view_grid_next();
        }
    }

    {
        plot::eps p("out/test_text.eps", 300, 200);
        p.font(24);
        p.view_inset(20);
        p.rmove_to(0.5,0.5).move(0,20).show("this^{shit^2} is a_1 |{test}", 0.5);
        p.rmove_to(0.5,0.5).move(0,-5).show("hello\nyou_{beta_{blop}}", 1);
    }

    {
        plot::eps p("out/test_plot.eps", 300, 200);
        p.font(24);
        p.view_inset(20);
        p.color(rgb(0.9,0.9,0.9));
        draw_axis_base(p);
        p.color(rgb(0,0,0));
        draw_axis(p, keywords(_btit("Mstar"), _ltit("SFR"), _ttit("Mstar"), _rtit("Lir")));
    }

    return 0;
}
