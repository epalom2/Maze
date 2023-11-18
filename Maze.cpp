#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>    // for std::stringstream
#include <stack>
#include <string>
#include <vector>

#include "catch.hpp"

// a simple type representing a position in the grid
struct grid_position
{
    int row;
    int column;

    friend bool operator==(grid_position const& lhs, grid_position const& rhs)
    {
        return lhs.row == rhs.row && lhs.column == rhs.column;
    }

    friend bool operator!=(grid_position const& lhs, grid_position const& rhs)
    {
        return !(lhs == rhs);
    }

    friend bool operator<(grid_position const& lhs, grid_position const& rhs)
    {
        return lhs.row < rhs.row ||
            (lhs.row == rhs.row && lhs.column < rhs.column);
    }

    // a > b --> !(b < a) || a == b
    friend bool operator>(grid_position const& lhs, grid_position const& rhs)
    {
        return !(rhs < lhs) || lhs == rhs;
    }

    // a >= b --> !(b > a)
    friend bool operator>=(grid_position const& lhs, grid_position const& rhs)
    {
        return !(lhs < rhs);
    }

    // a <= b ==> !(b < a)
    friend bool operator<=(grid_position const& lhs, grid_position const& rhs)
    {
        return !(rhs < lhs);
    }

    // input/output helper functions for grid_position
    friend std::istream& operator>>(std::istream& strm, grid_position& pos)
    {
        char c1, c2;
        strm >> c1 >> pos.row >> c2 >> pos.column;
        if (c1 != 'r' || c2 != 'c')
        {
            throw std::runtime_error("unrecognized format for grid_position");
        }
        return strm;
    }

    friend std::ostream& operator<<(
        std::ostream& strm, grid_position const& pos)
    {
        strm << 'r' << pos.row << 'c' << pos.column;
        return strm;
    }
};

// a simple grid type used to represent the maze
template <typename T>
class grid
{
public:
    grid(int size_x, int size_y, T init = T())
      : data(size_x * size_y, init)
      , columns(size_x)
      , rows(size_y)
    {
    }

    decltype(auto) operator[](grid_position const& pos)
    {
        return data[pos.column + pos.row * num_columns()];
    }
    decltype(auto) operator[](grid_position const& pos) const
    {
        return data[pos.column + pos.row * num_columns()];
    }

    int num_columns() const
    {
        return columns;
    }
    int num_rows() const
    {
        return rows;
    }

private:
    std::vector<T> data;
    int columns;
    int rows;
};

// read a maze description from a file
grid<bool> read_maze_from_file(std::string const& filename)
{
    // open file and verify that it was found
    std::ifstream strm(filename);
    if (!strm.is_open())
    {
        throw std::runtime_error("can't open file: " + filename);
    }

    // first read first line of file to extract the sizes of the maze
    int rows = 0, columns = 0;
    strm >> rows >> columns;
    if (rows == 0 || columns == 0)
    {
        throw std::runtime_error("unrecognized file format: " + filename);
    }

    // now read the maze data character by character
    grid<bool> maze(columns, rows);
    char c;

    std::noskipws(strm);    // prevent stream from skipping whitespace
    strm >> c;              // skip '\n'

    for (int y = 0; y < rows; ++y)
    {
        for (int x = 0; x < columns; ++x)
        {
            strm >> c;
            maze[grid_position{y, x}] = (c == ' ') ? true : false;
        }
        strm >> c;    // skip '\n'
    }

    return maze;
}

PROVIDED_TEST("verify reading maze file")
{
    grid<bool> g = read_maze_from_file("../mazes/maze11x11.txt");
    CHECK(g.num_columns() == 11);
    CHECK(g.num_rows() == 11);

    // test two arbitrary grid elements
    CHECK(!g[grid_position{0, 0}]);
    CHECK(g[grid_position{1, 2}]);
}

// read a solution for a maze from a file
std::vector<grid_position> read_solution_from_file(std::string const& filename)
{
    // open file and verify that it was found
    std::ifstream strm(filename);
    if (!strm.is_open())
    {
        throw std::runtime_error("can't open file: " + filename);
    }

    // first read first line of file to extract the sizes of the maze
    int positions = 0;
    strm >> positions;
    if (positions == 0)
    {
        throw std::runtime_error("unrecognized file format: " + filename);
    }

    std::vector<grid_position> path;
    path.reserve(positions);

    while (--positions >= 0)
    {
        grid_position pos;
        if (!(strm >> pos))
        {
            throw std::runtime_error("unrecognized file format: " + filename);
        }
        path.push_back(pos);
    }

    return path;
}

PROVIDED_TEST("verify reading solutions file")
{
    std::vector<grid_position> path =
        read_solution_from_file("../mazes/maze11x11.solution");

    std::vector<grid_position> expected = {grid_position{1, 1},
        grid_position{1, 2}, grid_position{1, 3}, grid_position{1, 4},
        grid_position{1, 5}, grid_position{1, 6}, grid_position{1, 7},
        grid_position{1, 8}, grid_position{1, 9}, grid_position{2, 9},
        grid_position{3, 9}, grid_position{4, 9}, grid_position{5, 9},
        grid_position{6, 9}, grid_position{7, 9}, grid_position{8, 9},
        grid_position{9, 9}};

    CHECK(path == expected);
}

// print given maze to stream
void print_maze_line_to_stream(
    std::ostream& strm, grid<bool> const& maze, int y)
{
    for (int x = 0; x < maze.num_columns(); ++x)
    {
        strm << (maze[grid_position{y, x}] ? " " : (x & 1 ? "-" : "+"));
    }
    strm << "\n";

    for (int x = 0; x < maze.num_columns(); ++x)
    {
        strm << (maze[grid_position{y + 1, x}] ? " " : "|");
    }
    strm << "\n";
}

void print_to_stream(std::ostream& strm, grid<bool> const& maze)
{
    for (int y = 0; y < maze.num_rows() - 1; y += 2)
    {
        print_maze_line_to_stream(strm, maze, y);
    }

    for (int x = 0; x < maze.num_columns(); ++x)
    {
        strm << (maze[grid_position{maze.num_rows() - 1, x}] ?
                " " :
                (x & 1 ? "-" : "+"));
    }
    strm << "\n";
}

char const* maze11x11 = R"(
+-+-+-+-+-+
|         |
+-+-+-+-+ +
|       | |
+ +-+-+ + +
| |   | | |
+ + + + + +
| | |   | |
+ +-+-+-+ +
|         |
+-+-+-+-+-+
)";

PROVIDED_TEST("verify maze printing")
{
    std::stringstream strm;
    strm << '\n';    // account for additional leading '\n' in maze11x11

    print_to_stream(strm, read_maze_from_file("../mazes/maze11x11.txt"));
    CHECK(maze11x11 == strm.str());
}

// Given a maze represented as a `grid<bool>` and a current `grid_position`
// cur, this function returns a `std::set` of all valid moves from `cur`.
std::set<grid_position> generate_valid_moves(grid<bool> const& maze, grid_position const& cur)
{
    std::set<grid_position> valid_moves;

    // Define the four cardinal directions (N, S, E, W)
    const std::vector<std::pair<int, int>> directions = {{-1, 0}, {1, 0}, {0, 1}, {0, -1}};

    // Check each cardinal direction from the current position
    for (const auto& dir : directions)
    {
        int new_row = cur.row + dir.first;
        int new_column = cur.column + dir.second;

        // Check if the new position is within the grid bounds
        if (new_row >= 0 && new_row < maze.num_rows() && new_column >= 0 && new_column < maze.num_columns())
        {
            
            // Check if the new position is an open corridor (not a wall)
            if (maze[grid_position{new_row, new_column}])
            {
                // Add the valid move to the set
                valid_moves.insert(grid_position{new_row, new_column});
            }
        }
    }

    return valid_moves;
}

PROVIDED_TEST("verifying generating valid moves from given position")
{
    grid<bool> maze = read_maze_from_file("../mazes/maze11x11.txt");
    std::set<grid_position> valid_moves =
        generate_valid_moves(maze, grid_position{1, 1});

    // we know there is just one possible move
    CHECK(valid_moves.size() == 1);

    // current position should not be part of valid moves
    CHECK(valid_moves.find(grid_position{1, 1}) == valid_moves.end());

    // we know that grid_position{1, 2} should be a valid move
    CHECK(valid_moves.find(grid_position{1, 2}) != valid_moves.end());
}

STUDENT_TEST("Generate valid moves for a different position in the maze")
{
   {//test 1
        grid<bool> maze = read_maze_from_file("../mazes/maze11x11.txt");
        std::set<grid_position> valid_moves = generate_valid_moves(maze, grid_position{1, 1});
        CHECK(valid_moves.size() == 1); 
   }
   
   {//test 2
        grid<bool> maze = read_maze_from_file("../mazes/maze11x11.txt");
        std::set<grid_position> valid_moves = generate_valid_moves(maze, grid_position{4, 5});
        CHECK(valid_moves.size() == 2); 
   }

    {//test 3
        grid<bool> maze = read_maze_from_file("../mazes/maze11x11.txt");
        std::set<grid_position> valid_moves = generate_valid_moves(maze, grid_position{9, 8});
        CHECK(valid_moves.size() == 2); 
   }

}


// Find the number of occurrences of a given grid_position contained in a
// given solution path
int count_known_positions(
    grid_position const& pos, std::vector<grid_position> const& p)
{
  // Sort the vector of grid positions
    std::vector<grid_position> sorted_positions = p;
    std::sort(sorted_positions.begin(), sorted_positions.end());

    // Use std::equal_range to find the range of occurrences of pos
    auto range = std::equal_range(sorted_positions.begin(), sorted_positions.end(), pos);

    // Calculate the number of occurrences
    int num_occurrences = std::distance(range.first, range.second);
    return num_occurrences;
}

PROVIDED_TEST("testing of a given position comprises a valid next step")
{
    std::vector<grid_position> path =
        read_solution_from_file("../mazes/maze11x11.solution");

    CHECK(1 == count_known_positions(grid_position{1, 1}, path));
    CHECK(0 == count_known_positions(grid_position{0, 0}, path));
}

STUDENT_TEST("More testing of a given position comprises a valid next step")
{

    std::vector<grid_position> path =
    read_solution_from_file("../mazes/maze11x11.solution"); 
    //test 1
    CHECK(1 == count_known_positions(grid_position{1, 2}, path));
    //test2
    CHECK(0 == count_known_positions(grid_position{4, 5}, path));
    //test3
    CHECK(0 == count_known_positions(grid_position{9, 8}, path));
}


// The function `validate_path` completes successfully (returns `true`) if
// all of the criteria for a valid solution are met. If it instead detects
// that the path violates one of the constraints as outlined in the
// assignment description, `validate_path`  should return `false`.
bool validate_path(
    grid<bool> const& maze, std::vector<grid_position> const& path)
{
    // Criterion 1: The path is not empty
    if (path.empty()) {
        return false;
    }

    // Criterion 2: The path starts at the entry (upper left corner) of the maze
    if (path.front() != grid_position{1, 1}) {
        return false;
    }

    // Criterion 3: The path ends at the exit (lower right corner) of the maze
    if (path.back() != grid_position{maze.num_rows() - 2, maze.num_columns() - 2}) {
        return false;
    }

    // Criterion 4: Each location in the path is a valid move from the previous path location
    for (size_t i = 1; i < path.size(); ++i) {
        grid_position prev_pos = path[i - 1];
        grid_position cur_pos = path[i];
        std::set<grid_position> valid_moves = generate_valid_moves(maze, prev_pos);
        if (valid_moves.find(cur_pos) == valid_moves.end()) {
            return false;
        }
    }

    // Criterion 5: The path must not contain a loop
    for (const auto& pos : path) {
        if (count_known_positions(pos, path) != 1) {
            return false;
        }
    }

    // All criteria are met, the path is valid
    return true;
}

PROVIDED_TEST("verifying validity of given solution path")
{
    grid<bool> maze = read_maze_from_file("../mazes/maze11x11.txt");
    std::vector<grid_position> path =
        read_solution_from_file("../mazes/maze11x11.solution");

    CHECK(validate_path(maze, path));
}

STUDENT_TEST("MORE verifying validity of given solution path")
{
  
   {//test 1
    grid<bool> maze = read_maze_from_file("../mazes/maze11x11.txt");
    std::vector<grid_position> path = read_solution_from_file("../mazes/maze11x11.solution");
    CHECK(validate_path(maze, path));
   }
   {
    //test 2
    //Invalid solution path - end position mismatch
    grid<bool> mazeA = read_maze_from_file("../mazes/maze11x11.txt");
    std::vector<grid_position> invalid_path = {grid_position{1, 1}, grid_position{3, 3}, grid_position{6, 6}};
    CHECK(!validate_path(mazeA, invalid_path));
   }
    {
    //test 3
    //Invalid solution path - contains a loop
    grid<bool> mazeB = read_maze_from_file("../mazes/maze11x11.txt");
    std::vector<grid_position> invalid_path = {grid_position{1, 1}, grid_position{2, 1}, grid_position{1, 1}};
    CHECK(!validate_path(mazeB, invalid_path));
    }
    {//test 4
        grid<bool> maze = read_maze_from_file("../mazes/maze7x9.txt");
        std::vector<grid_position> path = read_solution_from_file("../mazes/maze7x9.solution");
        CHECK(validate_path(maze, path)); 
    }

}



// Solve the given maze by finding the shortest solution path for it.
// Use the breadth-first search algorithm.

std::vector<grid_position> solve_maze_bfs(grid<bool> const& maze)
{
    std::deque<std::vector<grid_position>> paths;
    paths.push_back({grid_position{1, 1}}); // Start with a length-one path containing the entry position

    while (!paths.empty())
    {
        std::vector<grid_position> current_path = paths.front();
        paths.pop_front();

        // If the current path ends at the exit, it is the solution
        if (current_path.back() == grid_position{maze.num_rows() - 2, maze.num_columns() - 2})
        {
            return current_path;
        }

        // Get viable neighbors (valid moves not in the current path)
        grid_position last_position = current_path.back();
        std::set<grid_position> neighbors = generate_valid_moves(maze, last_position);

        for (const auto& neighbor : neighbors)
        {
            // Check if the neighbor is not in the current path
            if (count_known_positions(neighbor, current_path) == 0)
            {
                // Create a new path by extending the current path with the neighbor
                std::vector<grid_position> new_path = current_path;
                new_path.push_back(neighbor);
                paths.push_back(new_path);
            }
        }
    }

    // It should never reach this point because the problem assumes the maze has a solution
    throw std::runtime_error("No solution found for the maze.");
}



PROVIDED_TEST("verifying solving maze using breadth-first search")
{
    grid<bool> maze = read_maze_from_file("../mazes/maze11x11.txt");
    std::vector<grid_position> expected_path =
        read_solution_from_file("../mazes/maze11x11.solution");

    std::vector<grid_position> path = solve_maze_bfs(maze);

    CHECK(validate_path(maze, path));
    CHECK(expected_path == path);
}

STUDENT_TEST("More verifying solving maze using breadth-first search")
{
    {
    grid<bool> maze = read_maze_from_file("../mazes/maze23x23.txt");
    std::vector<grid_position> expected_path =
        read_solution_from_file("../mazes/maze23x23.solution");

    // Test case 1: Verify that the solution path obtained using BFS is valid
    std::vector<grid_position> bfs_path = solve_maze_bfs(maze);
    CHECK(validate_path(maze, bfs_path));

    // Test case 2: Verify that the solution path obtained using BFS matches the expected path
    CHECK(expected_path == bfs_path);
    }
    {//test 3
        grid<bool> maze = read_maze_from_file("../mazes/maze7x9.txt");
    std::vector<grid_position> expected_path =
        read_solution_from_file("../mazes/maze7x9.solution");

    std::vector<grid_position> path = solve_maze_bfs(maze);

    CHECK(validate_path(maze, path));
    CHECK(expected_path == path);
    }
}

// Solve the given maze by finding a solution path for it. Use the
// depth-first search algorithm.
std::vector<grid_position> solve_maze_dfs(grid<bool> const& maze)
{
    std::vector<grid_position> path;
    std::stack<std::vector<grid_position>> paths; // Stack to store paths

    // Start with a length-one path containing the entry position
    paths.push({grid_position{1, 1}});

    while (!paths.empty())
    {
        std::vector<grid_position> current_path = paths.top();
        paths.pop();

        // If the current path ends at the exit, it is the solution
        if (current_path.back() == grid_position{maze.num_rows() - 2, maze.num_columns() - 2})
        {
            return current_path;
        }

        // Get viable neighbors (valid moves not in the current path)
        grid_position last_position = current_path.back();
        std::set<grid_position> neighbors = generate_valid_moves(maze, last_position);

        for (const auto& neighbor : neighbors)
        {
            // Check if the neighbor is not in the current path
            if (count_known_positions(neighbor, current_path) == 0)
            {
                // Create a new path by extending the current path with the neighbor
                std::vector<grid_position> new_path = current_path;
                new_path.push_back(neighbor);
                paths.push(new_path);
            }
        }
    }

    // It should never reach this point because the problem assumes the maze has a solution
    throw std::runtime_error("No solution found for the maze.");
}


PROVIDED_TEST("verifying solving maze using depth-first search")
{
    grid<bool> maze = read_maze_from_file("../mazes/maze11x11.txt");
    std::vector<grid_position> expected_path =
        read_solution_from_file("../mazes/maze11x11.solution");

    std::vector<grid_position> path = solve_maze_dfs(maze);

    CHECK(validate_path(maze, path));
    CHECK(expected_path == path);
}

STUDENT_TEST("More verifying solving maze using depth-first search")
{
    {
    grid<bool> maze = read_maze_from_file("../mazes/maze23x23.txt");
    std::vector<grid_position> expected_path =
    read_solution_from_file("../mazes/maze23x23.solution");

    // Test case 1: Verify that the solution path obtained using DFS is valid
    std::vector<grid_position> dfs_path = solve_maze_dfs(maze);
    CHECK(validate_path(maze, dfs_path));

    // Test case 2: Verify that the solution path obtained using DFS matches the expected path
    CHECK(expected_path == dfs_path);
    }
    {
    //test 3
    grid<bool> maze = read_maze_from_file("../mazes/maze7x9.txt");
    std::vector<grid_position> expected_path =
    read_solution_from_file("../mazes/maze7x9.solution");
        
    std::vector<grid_position> path = solve_maze_dfs(maze);

    CHECK(validate_path(maze, path));
    CHECK(expected_path == path);
    }
}

void print_maze_line_with_solution_to_stream(
    std::ostream& strm, grid<bool> const& maze,
    std::vector<grid_position> const& solution, int y)
{
    for (int x = 0; x < maze.num_columns(); ++x)
    {
        grid_position current_position = {y, x};
        bool is_path = std::find(solution.begin(), solution.end(), current_position) != solution.end();
        strm << (maze[current_position] ? (is_path ? "*" : " ") : (x & 1 ? "-" : "+"));
    }
    strm << "\n";

    for (int x = 0; x < maze.num_columns(); ++x)
    {
        grid_position current_position = {y + 1, x};
        bool is_path = std::find(solution.begin(), solution.end(), current_position) != solution.end();
        strm << (maze[current_position] ? (is_path ? "*" : " ") : "|");
    }
    strm << "\n";
}

void print_maze_with_solution_to_stream(
    std::ostream& strm, grid<bool> const& maze,
    std::vector<grid_position> const& solution)
{
    strm << "\n";
    for (int y = 0; y < maze.num_rows() - 1; y += 2)
    {
        print_maze_line_with_solution_to_stream(strm, maze, solution, y);
    }

    for (int x = 0; x < maze.num_columns(); ++x)
    {
        grid_position current_position = {maze.num_rows() - 1, x};
        bool is_path = std::find(solution.begin(), solution.end(), current_position) != solution.end();
        strm << (maze[current_position] ? (is_path ? "*" : " ") : (x & 1 ? "-" : "+"));
    }
    strm << "\n";
}


STUDENT_TEST("Printing maze with solution path")
{
    // Define the maze and its solution path
    grid<bool> maze = read_maze_from_file("../mazes/maze11x11.txt");
    std::vector<grid_position> solution = read_solution_from_file("../mazes/maze11x11.solution");

    // Redirect standard output to a stringstream
    std::stringstream output;
    auto cout_buffer = std::cout.rdbuf(output.rdbuf());

    // Call the function to print maze with the solution path
    print_maze_with_solution_to_stream(output, maze, solution); // Use 'output' instead of 'std::cout'

    // Restore standard output
    std::cout.rdbuf(cout_buffer);

    // Define the expected output with the solution path overlaid
    char const* expected_output = R"(
+-+-+-+-+-+
|*********|
+-+-+-+-+*+
|       |*|
+ +-+-+ +*+
| |   | |*|
+ + + + +*+
| | |   |*|
+ +-+-+-+*+
|        *|
+-+-+-+-+-+
)";


// Compare the generated output with the expected output (ignoring leading/trailing whitespace)
CHECK((output.str()) == (expected_output));

}
