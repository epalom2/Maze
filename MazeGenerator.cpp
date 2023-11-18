#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

enum direction
{
    none,
    north = 1,
    east = 2,
    south = 4,
    west = 8
};

class maze
{
public:
    explicit maze(int size)
      : size(size)
      , pos_x(0)
      , pos_y(0)
      , data(size * size)
    {
    }

    void generate(std::string const& filename)
    {
        generate_maze();

        if (filename == "-")
        {
            write_to_file(std::cout);
        }
        else
        {
            std::ofstream strm(filename);
            write_to_file(strm);
        }
    }

private:
    void generate_maze()
    {
        while (true)
        {
            int d = get_direction();
            if (d < direction::north)
            {
                break;
            }

            switch (d)
            {
            case direction::north:
                cell(pos_x, pos_y) |= direction::north;
                cell(pos_x, --pos_y) = direction::south | direction::south << 4;
                break;
            case direction::east:
                cell(pos_x, pos_y) |= direction::east;
                cell(++pos_x, pos_y) = direction::west | direction::west << 4;
                break;
            case direction::south:
                cell(pos_x, pos_y) |= direction::south;
                cell(pos_x, ++pos_y) = direction::north | direction::north << 4;
                break;
            case direction::west:
                cell(pos_x, pos_y) |= direction::west;
                cell(--pos_x, pos_y) = direction::east | direction::east << 4;
            }
        }
    }

    void write_line(std::ostream& strm, int y) const
    {
        strm << "+";
        for (int x = 0; x < size; ++x)
        {
            strm << ((cell(x, y) & direction::north) ? " " : "+") << "+";
        }
        strm << "\n";

        strm << "|";
        for (int x = 0; x < size; ++x)
        {
            strm << ((cell(x, y) & direction::east) ? "  " : " |");
        }
        strm << "\n";
    }

    void write_to_file(std::ostream& strm) const
    {
        for (int y = 0; y != size; ++y)
        {
            write_line(strm, y);
        }

        strm << "+";
        for (int x = 0; x < size; ++x)
        {
            strm << ((cell(x, size - 1) & direction::south) ? " " : "-") << "+";
        }
        strm << "\n";
    }

    std::uint8_t& cell(int x, int y)
    {
        return data[x + size * y];
    }

    std::uint8_t cell(int x, int y) const
    {
        return data[x + size * y];
    }

    int get_direction()
    {
        int d = 1 << std::rand() % 4;
        while (true)
        {
            for (int x = 0; x < 4; x++)
            {
                if (test_direction(d))
                    return d;
                d <<= 1;
                if (d > direction::west)
                    d = direction::north;
            }

            d = (cell(pos_x, pos_y) & 0xf0) >> 4;
            if (d == direction::none)
                return d;

            switch (d)
            {
            case direction::north:
                --pos_y;
                break;
            case direction::east:
                ++pos_x;
                break;
            case direction::south:
                ++pos_y;
                break;
            case direction::west:
                --pos_x;
                break;
            }
            d = 1 << rand() % 4;
        }
        return direction::north;
    }

    bool test_direction(int d) const
    {
        switch (d)
        {
        case direction::north:
            return pos_y - 1 > -1 && !cell(pos_x, pos_y - 1);
        case direction::east:
            return pos_x + 1 < size && !cell(pos_x + 1, pos_y);
        case direction::south:
            return pos_y + 1 < size && !cell(pos_x, pos_y + 1);
        case direction::west:
            return pos_x - 1 > -1 && !cell(pos_x - 1, pos_y);
        }
        return false;
    }

private:
    int size;
    int pos_x, pos_y;
    std::vector<std::uint8_t> data;
};

bool odd(int n)
{
    return n & 1;
}

int main()
{
    std::srand(std::time(nullptr));

    std::cout
        << "Enter the maze size, an odd number bigger than 2 (0 to QUIT): ";

    int s = 0;
    std::cin >> s;
    if (!s)
    {
        return 0;    // quit
    }

    // make sure maze size is odd
    if (!odd(s))
    {
        ++s;
    }

    // make sure maze size is at least 3
    if (s < 3)
    {
        s = 3;
    }

    std::cout << "Enter the file name where to write the generated maze to (- "
                 "for std::cout): ";

    std::string filename("-");
    std::cin >> filename;
    if (filename.empty())
    {
        filename = "-";
    }

    maze m(s);
    m.generate(filename);

    std::cout << std::endl;
    return 0;
}
