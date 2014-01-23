#ifndef PLOT_HPP
#define PLOT_HPP

#include "color.hpp"
#include "math.hpp"
#include "string.hpp"
#include "print.hpp"

namespace plot {
    // Different stroke drawing styles for stroke endings
    enum class cap {
        def = 0,      // rough ending, exactly at the extremity of the line [default]
        round = 1,    // round ending, centered on the extremity of the line
        extended = 2  // same as 'def', but extended by 1/2 of the width of the stroke
    };

    // Different stroke drawing styles for stroke joins
    enum class join {
        def = 0,   // prolongate the lines until intersection [default]
        round = 1, // round join
        bevel = 2  // same as 'def', but truncated at line endings
    };

    // Helper class to write EPS plots
    class eps {
        friend class path;

        std::ofstream file_;
        uint_t font_size_ = 0;
        uint_t width_, height_;
        std::string default_font_ = "Helvetica";

        struct view_t {
            // TODO: switch to a matrix here
            // vec2d mat;
            float x0, y0, x1, y1;
            uint_t nx = 1, ny = 1, n = 0;
        };

        std::vector<view_t> views_;
        view_t* view_ = nullptr;


        void transform_abs_(float& x, float& y) {
            x += view_->x0;
            y += view_->y0;
        }

        void transform_rel_(float& x, float& y) {}

        void transform_scaled_abs_(float& x, float& y) {
            x *= (view_->x1 - view_->x0);
            y *= (view_->y1 - view_->y0);
            x += view_->x0;
            y += view_->y0;
        }

        void transform_scaled_rel_(float& x, float& y) {
            x *= (view_->x1 - view_->x0);
            y *= (view_->y1 - view_->y0);
        }

        void define_(const std::string& var, const std::string& value) {
            file_ << "/" << var << " (" << value << ") def\n";
        }

        void define_(const std::string& var, float value) {
            file_ << "/" << var << " " << value << " def\n";
        }

        struct fstring {
            std::string str;
            std::string font = "";
            int_t level = 0;
        };

        void parse_line_(const std::string& s, std::vector<fstring>& r, int_t level = 0) {
            fstring fs;
            fs.level = level;

            uint_t lpos = 0;
            uint_t opos = 0;
            uint_t pos = 0;
            while ((pos = s.find_first_of("^_|", opos)) != s.npos) {
                switch (s[pos]) {
                    case '|' : {
                        if (pos+1 != s.size() && s[pos+1] == '{') {
                            uint_t op = 1;
                            uint_t tp = pos+1;
                            while ((tp = s.find_first_of("{}", tp+1)) != s.npos) {
                                if (s[tp] == '{') ++op;
                                if (s[tp] == '}') --op;
                                if (op == 0) break;
                            }

                            if (op != 0) {
                                warning("parsing '", s,"':", pos+1, ": brace not closed");
                                opos = pos+1;
                                break;
                            }

                            fs.str += s.substr(lpos, pos - lpos);
                            r.push_back(fs);

                            // TODO: convert code to actual math character
                            fs.str = "<"+s.substr(pos+2, tp - pos-2)+">";
                            r.push_back(fs);
                            fs.str = "";

                            opos = tp+1;
                            lpos = tp+1;
                        } else {
                            fs.str += s.substr(lpos, pos - lpos);
                            lpos = pos+1;
                            opos = pos+2;
                        }
                        break;
                    }
                    case '^' :
                    case '_' : {
                        fs.str += s.substr(lpos, pos - lpos);
                        r.push_back(fs);
                        fs.str = "";

                        if (pos+1 == s.size()) {
                            warning("parsing '", s,"':", pos+1, ": trailing '", s[pos], "'");
                            opos = pos+1;
                            break;
                        }

                        std::string sub;
                        if (s[pos+1] == '{') {
                            uint_t op = 1;
                            uint_t tp = pos+1;
                            while ((tp = s.find_first_of("{}", tp+1)) != s.npos) {
                                if (s[tp] == '{') ++op;
                                if (s[tp] == '}') --op;
                                if (op == 0) break;
                            }

                            if (op != 0) {
                                warning("parsing '", s,"':", pos+1, ": brace not closed");
                                opos = pos+1;
                                break;
                            }

                            sub = s.substr(pos+2, tp - pos-2);

                            opos = tp+1;
                            lpos = tp+1;
                        } else {
                            sub = s[pos+1];
                            opos = pos+2;
                            lpos = pos+2;
                        }

                        parse_line_(sub, r, (s[pos] == '^' ? level+1 : level-1));
                        break;
                    }
                    default : {
                        opos = pos+1;
                        break;
                    }
                }
            }

            if (pos == s.npos) pos = s.size();
            fs.str += s.substr(lpos, pos - lpos);
            r.push_back(fs);
        }

        float level_scale_(int_t level) {
            return pow(0.5, abs(level));
        }

        float level_offset_(int_t level) {
            float f = 0.0;
            if (level > 0) {
                for (int_t i = 0; i < level; ++i) {
                    f += pow(0.7*0.5, i+1)/0.5;
                }
            } else {
                for (int_t i = 0; i < -level; ++i) {
                    f -= pow(0.4*0.5, i+1)/0.7;
                }
            }
            return f;
        }

        void show_line_(const std::string& text) {
            std::vector<fstring> fss;
            parse_line_(text, fss);
            int_t last_level = 0;
            uint_t size = font_size_;
            float height = char_height();
            
            for (auto& fs : fss) {
                if (fs.str.empty()) continue;

                if (fs.level != last_level) {
                    float doffset = level_offset_(fs.level) - level_offset_(last_level);
                    file_ << "0 " << height*doffset << " rmoveto\n";
                    font(size*level_scale_(fs.level));
                    last_level = fs.level;
                }

                file_ << "(" << fs.str << ") show\n";
            }

            if (last_level != 0) {
                float doffset = level_offset_(0) - level_offset_(last_level);
                file_ << "0 " << height*doffset << " rmoveto\n";
                font(size);
            }
        }
        
        void show_line_(const std::string& text, float align) {
            std::vector<fstring> fss;
            parse_line_(text, fss);
            int_t last_level = 0;
            uint_t size = font_size_;
            uint_t itxt = 0;

            for (auto& fs : fss) {
                if (fs.str.empty()) continue;

                if (fs.level != last_level) {
                    font(size*level_scale_(fs.level));
                    last_level = fs.level;
                }

                ++itxt;
                define_("text"+strn(itxt), fs.str);
                file_ << "text"+strn(itxt)+" stringwidth pop";
                if (itxt == 1) {
                    file_ << "\n";
                } else {
                    file_ << " add\n";
                }
            }

            if (last_level != 0) {
                font(size);
            }

            file_ << -clamp(align, 0, 1) << " mul 0 rmoveto\n";

            last_level = 0;
            size = font_size_;
            itxt = 0;
            float height = char_height();
            
            for (auto& fs : fss) {
                if (fs.str.empty()) continue;

                if (fs.level != last_level) {
                    float doffset = level_offset_(fs.level) - level_offset_(last_level);
                    file_ << "0 " << height*doffset << " rmoveto\n";
                    font(size*level_scale_(fs.level));
                    last_level = fs.level;
                }

                ++itxt;
                file_ << "text"+strn(itxt)+" show\n";
            }

            if (last_level != 0) {
                float doffset = level_offset_(0) - level_offset_(last_level);
                file_ << "0 " << height*doffset << " rmoveto\n";
                font(size);
            }
        }

    public :

    // -- Initialization and finalization
    // To create an EPS file, simply create a new object of the 'eps' class, give it a file name and
    // two integers for its dimensions (in points), then either let the object be destroyed or call
    // finish to finalize and save the file, once all the ploting is over.
    // Example: {
    //     plot::eps p("test.eps", 100, 100);
    //     // plot many things ...
    // }

        // Constructor: give file name and dimensions
        eps(const std::string& filename, uint_t width, uint_t height) :
            file_(filename), width_(width), height_(height) {

            view_reset();

            file_ << "%%!PS-Adobe-3.0 EPSF-3.0\n";
            file_ << "%%Creator: phypp\n";
            file_ << "%%DocumentData: Clean7Bit\n";
            file_ << "%%Origin: 0 0\n";
            file_ << "%%BoundingBox: 0 0 " << width_ << " " << height_ << "\n";
            file_ << "%%LanguageLevel: 2\n";
            file_ << "%%Pages: 1\n";
            file_ << "%%Page: 1 1\n\n";
        }

        // Destructor: closes the file properly (finish)
        ~eps() {
            finish();
        }

        // Close the file properly
        void finish() {
            file_ << "\n%%EOF\n";
            file_.close();
        }

    // -- Views into the file
    // View are a convenience feature that temporarily redefines the "apparent" origin and dimension
    // of the file, so that all subsequent plot calls will happen "as if" the origin was translated
    // and/or the dimensions of the file smaller. To do so, one has to define a "view box" inside
    // the file, using four coordinates x0, y0, x1 and y1. The class then takes care of transforming
    // any coordinate that is passed to 'move_to' or 'move' etc. By default, the view occupies the
    // whole file area (0, 0, with, height). This class handles a stack of views, so one can push a
    // new view *inside the current one*, and then come back to the previous one.

        // Return the width of the current view
        float width() const {
            return view_->x1 - view_->x0;
        }

        // Return the height of the current view
        float height() const {
            return view_->y1 - view_->y0;
        }

        // Return the horizontal origin of the current view
        float origin_x() const {
            return view_->x0;
        }

        // Return the vertical origin of the current view
        float origin_y() const {
            return view_->y0;
        }

        // Push an arbitrary view (x_left, y_bottom, x_right, y_top)
        eps& view(float x0, float y0, float x1, float y1) {
            view_t v;
            v.x0 = x0 + view_->x0; v.y0 = y0 + view_->y0;
            v.x1 = x1 + view_->x0; v.y1 = y1 + view_->y0;
            views_.push_back(v);
            view_ = &views_.back();
            return *this;
        }

        // Push an arbitrary scaled view (x_left, y_bottom, x_right, y_top)
        eps& rview(float x0, float y0, float x1, float y1) {
            view_t v;
            v.x0 = x0*(view_->x1 - view_->x0) + view_->x0;
            v.y0 = y0*(view_->y1 - view_->y0) + view_->y0;
            v.x1 = x1*(view_->x1 - view_->x0) + view_->x0;
            v.y1 = y1*(view_->y1 - view_->y0) + view_->y0;
            views_.push_back(v);
            view_ = &views_.back();
            return *this;
        }

        // Pop the last view from the stack
        eps& view_pop() {
            phypp_check(views_.size() != 1, "plot::eps: cannot pop the main view");
            views_.pop_back();
            view_ = &views_.back();
            return *this;
        }

        // Pop all views from the stack, and return to the default view (full file)
        eps& view_reset() {
            view_t v;
            v.x0 = 0.0f; v.y0 = 0.0f;
            v.x1 = width_; v.y1 = height_;
            views_.push_back(v);
            view_ = &views_.back();
            return *this;
        }

        // Create a path surrounding the current view
        eps& view_path() {
            start_path().move_to(0,0).rline_to(1,0).rline_to(1,1).rline_to(0,1).close_path();
            return *this;
        }

        // View the file as a grid of 'nx' horizontal regions and 'ny' vertical regions, and go the
        // 'n'th one:   #-------#
        //              | 2 | 3 |
        //              #---#---# 
        //              | 0 | 1 |
        //              #---#---# 
        eps& view_grid(int_t n, int_t nx, int_t ny) {
            float dx = width()/nx, dy = height()/ny;
            int_t x = n % ny, y = n / nx;

            view_t v;
            v.nx = nx; v.ny = ny; v.n = n % (nx*ny);
            v.x0 = x*dx + origin_x();     v.y0 = y*dy + origin_y();
            v.x1 = (x+1)*dx + origin_x(); v.y1 = (y+1)*dy + origin_y();
            views_.push_back(v);
            view_ = &views_.back();
            return *this;
        }

        // Shortcut for view_grid(0, nx, ny)
        eps& view_grid(int_t nx, int_t ny) {
            return view_grid(0, nx, ny);
        }

        // Shortcut for view_grid(n, previous_nx, previous_ny)
        // Note: this function does not add a new view on the stack, but replaces the last one
        eps& view_grid(int_t n) {
            view_t old = *view_;
            view_pop();
            return view_grid(n, old.nx, old.ny);
        }

        // Shortcut for view_grid(previous_n + 1, previous_nx, previous_ny)
        // Note: this function does not add a new view on the stack, but replaces the last one.
        // If using this function while being on the "last" element of the view grid, then the grid
        // view is completely popped, and the last view is restored.
        eps& view_grid_next() {
            view_t old = *view_;
            view_pop();
            if (old.n+1 != old.nx*old.ny) {
                return view_grid(old.n+1, old.nx, old.ny);
            } else {
                return *this;
            }
        }

        // Create a new view that is inset from the previous one (positive = inside).
        // 'ix0' : left inset, 'iy0' : bottom inset, 'ix1' : right inset, 'iy1' : top inset
        eps& view_inset(float ix0, float iy0, float ix1, float iy1) {
            return view(ix0, iy0, width() - ix1, height() - iy1);
        }

        // Create a new view that is inset from the previous one (positive = inside) (scaled
        // coordinates version).
        // 'ix0' : left inset, 'iy0' : bottom inset, 'ix1' : right inset, 'iy1' : top inset
        eps& rview_inset(float ix0, float iy0, float ix1, float iy1) {
            return rview(ix0, iy0, 1.0 - ix1, 1.0 - iy1);
        }

        // Shortcut for view_inset(i, i, i, i)
        eps& view_inset(float i) {
            return view_inset(i, i, i, i);
        }

        // Shortcut for view_rinset(i, i, i, i)
        eps& rview_inset(float i) {
            return rview_inset(i, i, i, i);
        }

    // -- Movement
    // The following functions allows one to move the current plot position inside the file,
    // either in an absolute sense ('move_to') or in a offset sense ('move'). One can also work in
    // scaled coordinates (0 to 1) using 'rmove_to' and 'rmove'. Note that these functions are
    // affected by the current view (see above).

        // Move to an absolute position in points
        eps& move_to(float x, float y) {
            transform_abs_(x, y);
            file_ << x << " " << y << " moveto\n";
            return *this;
        }

        // Move by a relative amount of points
        eps& move(float x, float y) {
            transform_rel_(x, y);
            file_ << x << " " << y << " rmoveto\n";
            return *this;
        }

        // Move to an absolute position in points
        eps& rmove_to(float x, float y) {
            transform_scaled_abs_(x, y);
            file_ << x << " " << y << " moveto\n";
            return *this;
        }

        // Move by a relative amount of points
        eps& rmove(float x, float y) {
            transform_scaled_rel_(x, y);
            file_ << x << " " << y << " rmoveto\n";
            return *this;
        }

        // Rotate the current coordinate system counterclockwise [degree]
        eps& rotate(float alpha) {
            file_ << alpha << " rotate\n";
            return *this;
        }

    // -- Graphic state
    // It is possible in EPS to save the current state (color, position, etc.) and restore it later.
    // This is done thanks to the two 'save' and 'restore' functions, whose name are quite self
    // explanatory. Note that is is only possible to save 32 different states at most at the same
    // time, and that restoring a state removes it from the saved state stack.

        // Save the current graphic state and put it on the state stack
        eps& save() {
            file_ << "gsave\n";
            return *this;
        }

        // Restore the most recently saved graphic state and pop it from the state stack
        eps& restore() {
            file_ << "grestore\n";
            return *this;
        }

    // -- Paths
    // Most of the drawing in EPS is done with paths. Paths are a defined as a set of points that
    // are used to draw a shape on the file. They are created point by point, can be disjoint or
    // fully closed. This class provides a proxy structure 'path' that exposes the path creation
    // functions 'line_to', 'line' and 'close_path', otherwise not available in the 'eps' class.
    // Also note that, as for 'move' and 'move_to', the 'line' and 'line_to' functions are affected
    // by the current view.
    // Example : {
    //     plot::eps p("test.eps", 100, 100);
    //     p.move_to(10,10);
    //     p.start_path().line(10,5).line(-10,5);
    //     // or:
    //     plot::eps::path(p).line(10,5).line(-10,5);
    // }

        // Path creation proxy
        struct path {
            eps& e;

            // Create a new path
            explicit path(eps& e_) : e(e_) {
                e.file_ << "newpath\n";
            }

            // Same as eps::move_to
            path& move_to(float x, float y) {
                e.move_to(x, y);
                return *this;
            }

            // Same as eps::move
            path& move(float x, float y) {
                e.move(x, y);
                return *this;
            }

            // Same as eps::rmove_to
            path& rmove_to(float x, float y) {
                e.rmove_to(x, y);
                return *this;
            }

            // Same as eps::rmove
            path& rmove(float x, float y) {
                e.rmove(x, y);
                return *this;
            }

            // Create a line from the current position to another absolute position
            path& line_to(float x, float y) {
                e.transform_abs_(x, y);
                e.file_ << x << " " << y << " lineto\n";
                return *this;
            }

            // Create a line from the current position to a relative offset
            path& line(float x, float y) {
                e.transform_rel_(x, y);
                e.file_ << x << " " << y << " rlineto\n";
                return *this;
            }

            // Create a line from the current position to another absolute position
            path& rline_to(float x, float y) {
                e.transform_scaled_abs_(x, y);
                e.file_ << x << " " << y << " lineto\n";
                return *this;
            }

            // Create a line from the current position to a relative offset
            path& rline(float x, float y) {
                e.transform_scaled_rel_(x, y);
                e.file_ << x << " " << y << " rlineto\n";
                return *this;
            }

            // Create a line from the current position to the first point of the path
            path& close_path() {
                e.file_ << "closepath\n";
                return *this;
            }

            // Convenience function to return to the 'eps' file, calling this function is not needed
            eps& end_path() {
                return e;
            }
        };

        // Create a new path
        path start_path() {
            return path(*this);
        }

    // -- Path drawing
    // Once a path has been created, it will not display anything in the file. One has to call some
    // drawing functions such as the ones below in order for the path to be visible. In particular,
    // it is possible to draw the path as a stroke ('stroke') or to fill it with some color
    // ('fill'). Note however that such action destroys the path, so if one needs to plot both a
    // stroke and fill the path, one has to save the path after it was created, call 'stroke', 
    // restore the path, then call 'fill'.

        // Set the current drawing color in RGB format
        eps& color(const rgb& c) {
            file_ << c.r << " " << c.g << " " << c.b << " setrgbcolor\n";
            return *this;
        }

        // Set the width to use when drawing strokes (in points)
        eps& line_width(float w) {
            file_ << w << " setlinewidth\n";
            return *this;
        }

        // Set the dash pattern
        // The pattern specifies an alternance of plain / transparent dashes, so that [a b c d] will
        // produce a pattern with 'a' points of plain stroke, 'b' points of empty space, 'c' point
        // of plain stroke again, then 'd' points of empty space once more, etc. Specifying an empty
        // pattern goes back to plain stroke drawing.
        eps& line_dashed(const vec1f& pattern, float offset = 0.0) {
            file_ << "[";
            for (uint_t i = 0; i < pattern.size(); ++i) {
                if (i != 0) file_ << " ";
                file_ << pattern[i];
            }
            file_ << "] " << offset << " setdash\n";
            return *this;
        }

        // Remove any dash pattern and go back to plain stroke drawing
        eps& line_plain() {
            return line_dashed({}, 0.0f);
        }

        // Set the stroke end style, see plot::cap
        eps& line_cap(cap c) {
            file_ << int(c) << " setlinecap\n";
            return *this;
        }

        // Set the stroke join style, see plot::join
        eps& line_join(join j) {
            file_ << int(j) << " setlinejoin\n";
            return *this;
        }

        // Stroke the current path and destroy it
        eps& stroke() {
            file_ << "stroke\n";
            return *this;
        }

        // Fill the current path and destroy it
        eps& fill() {
            file_ << "fill\n";
            return *this;
        }

    // -- Text drawing
    // Drawing text in EPS first require a font. Once one has been chosen (using 'font'), one can
    // display any amount of text using the 'show' function. Note that EPS does not natively support
    // line jumps, but this class can handle it by looking for manual line jumps '\n' in the string,
    // and manually going to the next line to continue drawing. To do so, it must use one save slot.

        // Set the current font for drawing text and the fond size (in point)
        eps& font(const std::string& name, uint_t size) {
            if (default_font_ == name && font_size_ == size) return *this;
            default_font_ = name;
            font_size_ = size;
            file_ << "/" << name << " findfont " << size << " scalefont setfont\n";
            return *this;
        }

        // Use the previous font with a modified size (in point)
        eps& font(uint_t size) {
            if (font_size_ == size) return *this;
            font_size_ = size;
            file_ << "/" << default_font_ << " findfont " << size << " scalefont setfont\n";
            return *this;
        }

        // Display some text at the current position (left aligned)
        eps& show(const std::string& text) {
            vec1s lines = split(text, "\n");
            if (lines.size() == 1) {
                show_line_(text);
            } else {
                save();
                for (uint_t i = 0; i < lines.size(); ++i) {
                    if (i != 0) {
                        restore();
                        file_ << "0 " << -float(i*font_size_) << " rmoveto\n";
                    }
                    show_line_(text);
                }
            }
            
            return *this;
        }

        // Display some text at the current position (aligned : 0.0 = leftmost, 1.0 = rightmost)
        eps& show(const std::string& text, float align) {
            if (align == 0.0f) return show(text);

            vec1s lines = split(text, "\n");
            if (lines.size() == 1) {
                show_line_(text, align);
            } else {
                save();
                for (uint_t i = 0; i < lines.size(); ++i) {
                    if (i != 0) {
                        restore();
                        file_ << "0 " << -float(i*font_size_) << " rmoveto\n";
                    }
                    show_line_(lines[i], align);
                }
            }
            
            return *this;
        }

        // Return the height of a character using the current font
        float char_height() const {
            return 0.7*font_size_;
        }

    // -- Clipping
    // It is possible to use a path to select a particular portion of the file, outside of which
    // nothing will be displayed by future drawing calls. This is particularly usefull for defining
    // plotting areas, where the data cannot be displayed outside of the plot's boundary. By
    // default, this clipping region is set to the whole file, so it is virtually inactive.

        // Destroy the current path and use it as a clipping region
        eps& clip() {
            file_ << "clip\n";
            return *this;
        }

        // Remove the clipping region
        eps& no_clip() {
            view_t* ov = view_;
            view_path().clip();
            view_ = ov;
            return *this;
        }
    };

    void draw_axis_base(eps& e) {
        e.view_path();
        e.stroke();
    }

    void draw_axis(eps& e, declare_keywords(_ltit(""), _btit(""), _rtit(""), _ttit(""))) {
        bool inset = false;
        float ix0 = 0.0f, iy0 = 0.0f, ix1 = 0.0f, iy1 = 0.0f;

        if (keyword_set(_btit)) {
            e.rmove_to(0.5,0);
            e.show(get_keyword(_btit), 0.5);
            iy0 = e.char_height()*1.3;
            inset = true;
        }

        if (keyword_set(_ttit)) {
            e.rmove_to(0.5,1).move(0,-e.char_height());
            e.show(get_keyword(_ttit), 0.5);
            iy1 = e.char_height()*1.3;
            inset = true;
        }

        if (keyword_set(_ltit)) {
            e.rmove_to(0,0.5).move(e.char_height(),0);
            e.rotate(90).show(get_keyword(_ltit), 0.5).rotate(-90);
            ix0 = e.char_height()*1.3;
            inset = true;
        }

        if (keyword_set(_rtit)) {
            e.rmove_to(1,0.5);
            e.rotate(90).show(get_keyword(_rtit), 0.5).rotate(-90);
            ix1 = e.char_height()*1.3;
            inset = true;
        }

        if (inset) {
            e.view_inset(ix0, iy0, ix1, iy1);
        }

        draw_axis_base(e);
        
        if (inset) {
            e.view_pop();
        }
    }
}

#endif
