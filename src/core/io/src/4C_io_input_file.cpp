// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_io_input_file.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_utils_string.hpp"

#include <ryml.hpp>
#include <ryml_std.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  namespace
  {

    /**
     * Sections that contain at least this number of entries are considered huge and are only
     * available on rank 0.
     */
    constexpr std::size_t huge_section_threshold = 10'000;

    //! The different ways we want to handle sections in the input file.
    enum class SectionType
    {
      //! A section that is read directly.
      normal,
      //! A section that mentions other files that are included and need to be read.
      include,
    };


    std::filesystem::path get_include_path(
        const std::string& include_line, const std::filesystem::path& current_file)
    {
      // Interpret the path as relative to the currently read file, if is not absolute
      std::filesystem::path included_file(include_line);
      if (!included_file.is_absolute())
      {
        included_file = current_file.parent_path() / included_file;
      }
      FOUR_C_ASSERT_ALWAYS(
          std::filesystem::status(included_file).type() == std::filesystem::file_type::regular,
          "Included file '%s' is not a regular file. Does the file exist?", included_file.c_str());
      return included_file;
    }

    void join_lines(std::list<std::string>& list_of_lines, std::vector<char>& raw_content,
        std::vector<std::string_view>& lines)
    {
      FOUR_C_ASSERT(raw_content.empty() && lines.empty(),
          "Implementation error: raw_content and lines must be empty.");

      // Sum up the length of all lines to reserve the memory for the raw content.
      const std::size_t raw_content_size =
          std::accumulate(list_of_lines.begin(), list_of_lines.end(), std::size_t{0},
              [](std::size_t sum, const auto& line) { return sum + line.size(); });

      raw_content.reserve(raw_content_size);
      lines.reserve(list_of_lines.size());

      for (const auto& line : list_of_lines)
      {
        FOUR_C_ASSERT(raw_content.data(), "Implementation error: raw_content must be allocated.");
        const auto* start_of_line = raw_content.data() + raw_content.size();
        raw_content.insert(raw_content.end(), line.begin(), line.end());
        lines.emplace_back(start_of_line, line.size());
      }
    }

    std::vector<std::filesystem::path> read_dat_content(const std::filesystem::path& file_path,
        std::unordered_map<std::string, InputFile::SectionContent>& content_by_section)
    {
      const auto name_of_section = [](const std::string& section_header)
      {
        auto pos = section_header.rfind("--");
        if (pos == std::string::npos) return std::string{};
        return Core::Utils::trim(section_header.substr(pos + 2));
      };

      std::ifstream file(file_path);
      if (not file) FOUR_C_THROW("Unable to open file: %s", file_path.c_str());

      // Tracking variables while walking through the file
      std::vector<std::filesystem::path> included_files;
      SectionType current_section_type = SectionType::normal;
      std::list<std::string> list_of_lines;
      InputFile::SectionContent* current_section_content = nullptr;
      std::string line;

      // Loop over all input lines. This reads the actual file contents and determines whether a
      // line is to be read immediately or should be excluded because it is in one of the excluded
      // sections.
      while (getline(file, line))
      {
        // In case we are reading an include section, a comment needs to be preceded by
        // whitespace. Otherwise, we would treat double slashes as comments, although they are
        // part of the file path.
        if (current_section_type == SectionType::include)
        {
          // Take care to remove comments only if they are preceded by whitespace.
          line = Core::Utils::strip_comment(line, " //");
          if (line.empty()) continue;

          // Additionally check if the first token is a comment to handle the case where the
          // comment starts at the beginning of the line.
          if (line.starts_with("//")) continue;
        }
        // Remove comments, trailing and leading whitespaces, compact internal whitespaces
        else
        {
          line = Core::Utils::strip_comment(line);
        }

        // line is now empty
        if (line.size() == 0) continue;

        // This line starts a new section
        if (line.starts_with("--"))
        {
          // Finish the current section.
          if (current_section_content)
          {
            join_lines(list_of_lines, current_section_content->raw_content,
                current_section_content->lines);
          }
          list_of_lines.clear();

          const auto name = name_of_section(line);
          FOUR_C_ASSERT(name.size() > 0, "Section name must not be empty.");

          // Determine what kind of new section we started.
          if (line.rfind("--INCLUDES") != std::string::npos)
          {
            current_section_type = SectionType::include;
            current_section_content = nullptr;
          }
          else
          {
            current_section_type = SectionType::normal;
            FOUR_C_ASSERT_ALWAYS(content_by_section.find(name) == content_by_section.end(),
                "Section '%s' is defined again in file '%s'.", name.c_str(), file_path.c_str());

            content_by_section[name] = {};
            current_section_content = &content_by_section[name];
            current_section_content->file = file_path;
          }
        }
        // The line is part of a section.
        else
        {
          switch (current_section_type)
          {
            case SectionType::normal:
            {
              list_of_lines.emplace_back(line);
              break;
            }
            case SectionType::include:
            {
              if (!line.starts_with("--"))
              {
                included_files.emplace_back(get_include_path(line, file_path));
              }
              break;
            }
          }
        }
      }

      // Finish the current section.
      if (current_section_content)
      {
        join_lines(
            list_of_lines, current_section_content->raw_content, current_section_content->lines);
      }

      return included_files;
    }

    std::vector<std::filesystem::path> read_yaml_content(const std::filesystem::path& file_path,
        std::unordered_map<std::string, InputFile::SectionContent>& content_by_section)
    {
      std::vector<std::filesystem::path> included_files;

      // In this first iteration of the YAML support, we map the constructs from a YAML file back
      // to constructs in a dat file. This means that top-level sections are pre-fixed with "--" and
      // the key-value pairs are mapped to "key = value" lines.
      //

      // Read the whole file into a string and parse it with ryml.
      std::ifstream file(file_path);
      std::ostringstream ss;
      ss << file.rdbuf();
      const std::string file_content = ss.str();
      ryml::Tree tree = ryml::parse_in_arena(ryml::to_csubstr(file_content));

      const auto to_string = [](const ryml::csubstr str) -> std::string
      { return std::string(str.data(), str.size()); };

      for (ryml::ConstNodeRef node : tree.crootref())
      {
        const auto& section_name = node.key();

        // If this is the special section "INCLUDES", we need to handle it differently.
        if (section_name == "INCLUDES")
        {
          if (node.has_val())
          {
            included_files.emplace_back(get_include_path(to_string(node.val()), file_path));
          }
          else if (node.is_seq())
          {
            for (const auto& include_node : node)
            {
              included_files.emplace_back(
                  get_include_path(to_string(include_node.val()), file_path));
            }
          }
          else
            FOUR_C_THROW("INCLUDES section must contain a single file or a sequence.");

          continue;
        }


        FOUR_C_ASSERT_ALWAYS(
            content_by_section.find(to_string(section_name)) == content_by_section.end(),
            "Section '%s' is defined again in file '%s'.", to_string(section_name).c_str(),
            file_path.c_str());

        auto& current_content = content_by_section[to_string(section_name)];
        current_content.file = file_path;
        std::list<std::string> list_of_lines;

        const auto read_flat_sequence = [&](const ryml::ConstNodeRef& node)
        {
          for (const auto& entry : node)
          {
            FOUR_C_ASSERT_ALWAYS(entry.has_val(),
                "While reading section '%s': "
                "only scalar entries are supported in sequences.",
                to_string(section_name).c_str());
            list_of_lines.emplace_back(to_string(entry.val()));
          }
        };

        const auto read_map = [&](const ryml::ConstNodeRef& node)
        {
          for (const auto& entry : node)
          {
            FOUR_C_ASSERT_ALWAYS(entry.has_val(),
                "While reading section '%s': "
                "only scalar key-value pairs are supported in maps.",
                to_string(section_name).c_str());
            list_of_lines.emplace_back(to_string(entry.key()) + " = " + to_string(entry.val()));
          }
        };

        if (node.is_map())
        {
          read_map(node);
        }
        else if (node.is_seq())
        {
          read_flat_sequence(node);
        }
        else
        {
          FOUR_C_THROW("Entries in section %s must either form a map or a sequence",
              to_string(section_name).c_str());
        }

        // Finish the current section by condensing the lines into the content.
        join_lines(list_of_lines, current_content.raw_content, current_content.lines);
      }

      return included_files;
    }
  }  // namespace

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  InputFile::InputFile(std::string filename, MPI_Comm comm) : comm_(comm)
  {
    read_generic(filename);
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  std::filesystem::path InputFile::file_for_section(const std::string& section_name) const
  {
    auto it = content_by_section_.find(section_name);
    if (it == content_by_section_.end()) return std::filesystem::path{};
    return std::filesystem::path{it->second.file};
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  bool InputFile::has_section(const std::string& section_name) const
  {
    const bool known_somewhere = Core::Communication::all_reduce<bool>(
        content_by_section_.contains(section_name),
        [](const bool& r, const bool& in) { return r || in; }, comm_);
    return known_somewhere;
  }



  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void InputFile::read_generic(const std::filesystem::path& top_level_file)
  {
    if (Core::Communication::my_mpi_rank(comm_) == 0)
    {
      // Start by "including" the top-level file.
      std::list<std::filesystem::path> included_files{top_level_file};

      // We use a hand-rolled loop here because the list keeps growing; thus we need to
      // continuously re-evaluate where the end of the list is.
      for (auto file_it = included_files.begin(); file_it != included_files.end(); ++file_it)
      {
        std::vector<std::filesystem::path> new_include_files = std::invoke(
            [&]
            {
              const auto file_extension = file_it->extension().string();
              // Note that json is valid yaml and we can read it with the yaml parser.
              if (file_extension == ".yaml" || file_extension == ".yml" ||
                  file_extension == ".json")
              {
                return read_yaml_content(*file_it, content_by_section_);
              }
              else
              {
                return read_dat_content(*file_it, content_by_section_);
              }
            });

        // Check that the file is not included twice
        for (const auto& file : new_include_files)
        {
          if (std::ranges::find(included_files, file) != included_files.end())
          {
            FOUR_C_THROW(
                "File '%s' was already included before.\n Cycles are not allowed.", file.c_str());
          }
          else
          {
            included_files.emplace_back(file);
          }
        }
      }
    }

    if (Core::Communication::my_mpi_rank(comm_) == 0)
    {
      // Temporarily move the sections that are not huge into a separate map.
      std::unordered_map<std::string, SectionContent> non_huge_sections;

      for (auto&& [section_name, content] : content_by_section_)
      {
        if (content.lines.size() < huge_section_threshold)
        {
          non_huge_sections[section_name] = std::move(content);
        }
      }

      Core::Communication::broadcast(non_huge_sections, 0, comm_);

      // Move the non-huge sections back into the main map.
      for (auto&& [section_name, content] : non_huge_sections)
      {
        content_by_section_[section_name] = std::move(content);
      }
    }
    else
    {
      // Other ranks receive the non-huge sections.
      Core::Communication::broadcast(content_by_section_, 0, comm_);
    }

    // the following section names are always regarded as valid
    record_section_used("TITLE");
    record_section_used("FUNCT1");
    record_section_used("FUNCT2");
    record_section_used("FUNCT3");
    record_section_used("FUNCT4");
    record_section_used("FUNCT5");
    record_section_used("FUNCT6");
    record_section_used("FUNCT7");
    record_section_used("FUNCT8");
    record_section_used("FUNCT9");
    record_section_used("FUNCT10");
    record_section_used("FUNCT11");
    record_section_used("FUNCT12");
    record_section_used("FUNCT13");
    record_section_used("FUNCT14");
    record_section_used("FUNCT15");
    record_section_used("FUNCT16");
    record_section_used("FUNCT17");
    record_section_used("FUNCT18");
    record_section_used("FUNCT19");
    record_section_used("FUNCT20");
  }



  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  bool InputFile::print_unknown_sections(std::ostream& out) const
  {
    using MapType = decltype(knownsections_);
    const MapType merged_map = Core::Communication::all_reduce<MapType>(
        knownsections_,
        [](const MapType& r, const MapType& in)
        {
          MapType result = r;
          for (const auto& [key, value] : in)
          {
            result[key] |= value;
          }
          return result;
        },
        comm_);
    const bool printout = std::any_of(
        merged_map.begin(), merged_map.end(), [](const auto& kv) { return !kv.second; });

    // now it's time to create noise on the screen
    if (printout and (Core::Communication::my_mpi_rank(get_comm()) == 0))
    {
      out << "\nERROR!"
          << "\n--------"
          << "\nThe following input file sections remained unused (obsolete or typo?):\n";
      for (const auto& [section_name, known] : knownsections_)
      {
        if (!known) out << section_name << '\n';
      }
      out << std::endl;
    }

    return printout;
  }


  void InputFile::record_section_used(const std::string& section_name)
  {
    knownsections_[section_name] = true;
  }


  void InputFile::SectionContent::pack(Communication::PackBuffer& data) const
  {
    Core::Communication::add_to_pack(data, file);
    Core::Communication::add_to_pack(data, raw_content);

    // String_views are not packable, so we store offsets.
    std::vector<std::size_t> offsets;
    offsets.reserve(lines.size());
    for (const auto& line : lines)
    {
      FOUR_C_ASSERT(line.data() >= raw_content.data() &&
                        line.data() <= raw_content.data() + raw_content.size(),
          "Line data out of bounds.");
      const std::size_t offset = (line.data() - raw_content.data()) / sizeof(char);
      FOUR_C_ASSERT(offset <= raw_content.size(), "Offset out of bounds.");
      offsets.push_back(offset);
    }
    FOUR_C_ASSERT(offsets.empty() || offsets.back() + lines.back().size() == raw_content.size(),
        "Offset out of bounds.");
    Core::Communication::add_to_pack(data, offsets);
  }


  void InputFile::SectionContent::unpack(Communication::UnpackBuffer& buffer)
  {
    Core::Communication::extract_from_pack(buffer, file);
    Core::Communication::extract_from_pack(buffer, raw_content);

    std::vector<std::size_t> offsets;
    Core::Communication::extract_from_pack(buffer, offsets);
    lines.clear();
    for (std::size_t i = 0; i < offsets.size(); ++i)
    {
      const char* start = raw_content.data() + offsets[i];
      const std::size_t length = (i + 1 < offsets.size() ? (offsets[i + 1] - offsets[i])
                                                         : (raw_content.size()) - offsets[i]);
      FOUR_C_ASSERT(
          start >= raw_content.data() && start + length <= raw_content.data() + raw_content.size(),
          "Line data out of bounds.");
      lines.emplace_back(start, length);
    }
  }


  std::pair<std::string, std::string> read_key_value(const std::string& line)
  {
    std::string::size_type separator_index = line.find('=');
    // The equals sign is only treated as a separator when surrounded by whitespace.
    if (separator_index != std::string::npos &&
        !(std::isspace(line[separator_index - 1]) && std::isspace(line[separator_index + 1])))
      separator_index = std::string::npos;

    // In case we didn't find an "=" separator, look for a space instead
    if (separator_index == std::string::npos)
    {
      separator_index = line.find(' ');

      if (separator_index == std::string::npos)
        FOUR_C_THROW("Line '%s' with just one word in parameter section", line.c_str());
    }

    std::string key = Core::Utils::trim(line.substr(0, separator_index));
    std::string value = Core::Utils::trim(line.substr(separator_index + 1));

    if (key.empty()) FOUR_C_THROW("Cannot get key from line '%s'", line.c_str());
    if (value.empty()) FOUR_C_THROW("Cannot get value from line '%s'", line.c_str());

    return {std::move(key), std::move(value)};
  }

}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE
