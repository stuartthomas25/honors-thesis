
#include <ostream>
#include <string>

class progress_bar
{
    static const auto overhead = sizeof " [100%]";

    ostream& os;
    const size_t bar_width;
    string message;
    const string full_bar;

 public:
    progress_bar(ostream& os, size_t line_width,
                 string message_, const char symbol = '.')
        : os{os},
          bar_width{line_width - overhead},
          message{move(message_)},

          full_bar{string(bar_width, symbol) + string(bar_width, ' ')}
    {
        if (message.size()+1 >= bar_width || message.find('\n') != message.npos) {
            os << message << '\n';
            message.clear();
        } else {
            message += ' ';
        }
        write(0.0);
    }

    // not copyable
    progress_bar(const progress_bar&) = delete;
    progress_bar& operator=(const progress_bar&) = delete;

    ~progress_bar()
    {
        write(1.0);
        //os << '\n';
    }

    void write(double fraction);
};

void progress_bar::write(double fraction)
{
    // clamp fraction to valid range [0,1]
    if (fraction < 0)
        fraction = 0;
    else if (fraction > 1)
        fraction = 1;

    auto width = bar_width - message.size();
    auto offset = bar_width - static_cast<unsigned>(width * fraction);

    os << '\r' << message;
    os.write(full_bar.data() + offset, width);
    os << " [" << setw(3) << static_cast<int>(100*fraction) << "%] " << flush;
} 

