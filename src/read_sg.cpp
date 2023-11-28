struct ReadInfo {
    int read_id;
    int loc;
    int read_length;
};

bool StrGraph::readSGFile(std::string& filename) {
    clock_t t = clock();
    std::ifstream fin(filename);

    if (!fin.is_open()) {
        std::cout << "Error: " << filename << " does not exist!\n";
        return false;
    }

    std::string line;
    boost::char_separator<char> sep(",");

    std::vector<ReadInfo> read_info;
    read_info.reserve(10000); // Reserve an initial capacity to reduce reallocations

    while (std::getline(fin, line)) {
        if (line[0] == '>') {
            line.erase(0, 1);
            boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
            size_t colon = line.find(':');

            if (colon == std::string::npos) {
                processNodeHeader(tokens, read_info);
            } else {
                processEdgeHeader(tokens, read_info);
            }
        }
    }

    p_order = p_seq.size();
    ff_stringGraph = true;
    fin.clear();
    fin.close();
    t = clock() - t;
    std::cout << "Load fastq takes " << static_cast<double>(t) / CLOCKS_PER_SEC << " s\n";

    return true;
}

void StrGraph::processNodeHeader(const boost::tokenizer<boost::char_separator<char>>& tokens, std::vector<ReadInfo>& read_info) {
    auto it = tokens.begin();
    int node_id = std::stoi(*it);
    p_all_nodes.push_back(node_id);
    ++it;
    int read_id = -1;
    IntegerVector read_on_this_node, loc_on_this_node;

    for (int i = 0; it != tokens.end(); ++it, ++i) {
        if (i % 3 == 0) {
            read_on_this_node.push_back(std::stoi(*it));
        } else if (i % 3 == 1) {
            loc_on_this_node.push_back(std::stoi(*it));
        } else {
            int read_length = std::stoi(*it);
            p_read_length[read_on_this_node.back()] = read_length;
            read_info.push_back({read_on_this_node.back(), loc_on_this_node.back(), read_length});
        }
    }

    p_read_on_node.push_back(read_on_this_node);
    p_pos_on_node.push_back(loc_on_this_node);
    BoostSTRVertex v_source;

    if (vertex.count(node_id) == 0) {
        STRVertexType node(node_id);
        v_source = boost::add_vertex(node, *p_graph_);
        vertex[node_id] = v_source;
    } else {
        v_source = vertex[node_id];
    }
}

void StrGraph::processEdgeHeader(const boost::tokenizer<boost::char_separator<char>>& tokens, std::vector<ReadInfo>& read_info) {
    auto it = tokens.begin();
    int node_id = std::stoi(*it);
    p_all_nodes.push_back(node_id);
    ++it;
    IntegerVector read_on_this_node, loc_on_this_node;

    for (int i = 0; it != tokens.end(); ++it, ++i) {
        if (i % 3 == 0) {
            read_on_this_node.push_back(std::stoi(*it));
        } else if (i % 3 == 1) {
            loc_on_this_node.push_back(std::stoi(*it));
        } else {
            int read_length = std::stoi(*it);
            p_read_length[read_on_this_node.back()] = read_length;
            read_info.push_back({read_on_this_node.back(), loc_on_this_node.back(), read_length});
        }
    }

    p_read_on_node.push_back(read_on_this_node);
    p_pos_on_node.push_back(loc_on_this_node);
    BoostSTRVertex v_source, v_target;

    if (vertex.count(node_id) == 0) {
        STRVertexType node(node_id);
        v_source = boost::add_vertex(node, *p_graph_);
        vertex[node_id] = v_source;
    } else {
        v_source = vertex[node_id];
    }

    // Parse the rest of the header and create edges
    processEdgeRest(tokens.begin(), tokens.end(), v_source, read_info);
}

void StrGraph::processEdgeRest(boost::tokenizer<boost::char_separator<char>>::iterator it, boost::tokenizer<boost::char_separator<char>>::iterator end, const BoostSTRVertex& v_source, std::vector<ReadInfo>& read_info) {
    // Parse and create edges based on the remaining tokens
    for (; it != end; ++it) {
        int target_id = std::stoi(*it);
        BoostSTRVertex v_target;

        if (vertex.count(target_id) > 0) {
            v_target = vertex[target_id];
        } else {
            STRVertexType node(target_id);
            v_target = boost::add_vertex(node, *p_graph_);
            vertex[target_id] = v_target;
        }

        std::pair<BoostSTREdge, bool> e_search = boost::add_edge(v_source, v_target, *p_graph_);

        if (!e_search.second) {
            std::cout << "Error: StringGraph::LoadGraph: Failed to add edges between vertices!\n";
        }
    }
}
