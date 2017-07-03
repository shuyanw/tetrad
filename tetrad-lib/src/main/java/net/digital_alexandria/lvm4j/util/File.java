/**
 * lvm4j: a Java implementation of various latent variable models.
 *
 * Copyright (C) 2015 - 2016 Simon Dirmeier
 *
 * This file is part of lvm4j.
 * <p>
 * lvm4j is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * lvm4j is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with lvm4j.  If not, see <http://www.gnu.org/licenses/>.
 */

package net.digital_alexandria.lvm4j.util;

import net.digital_alexandria.lvm4j.datastructures.Pair;
import net.digital_alexandria.lvm4j.datastructures.Triple;
import net.digital_alexandria.lvm4j.edges.WeightedArc;
import net.digital_alexandria.lvm4j.enums.ExitCode;
import net.digital_alexandria.lvm4j.markovmodel.HMM;
import net.digital_alexandria.lvm4j.markovmodel.HMMNode;
import net.digital_alexandria.lvm4j.markovmodel.HMMParams;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * @author Simon Dirmeier {@literal s@simon-dirmeier.net}
 */
public final class File
{

    private final static Logger _LOGGER = LoggerFactory.getLogger(File.class);

    private static final String xmlDefinition = "<markovmodel>\n" +
                                                 "\t<meta>\n" +
                                                 "\t\t<states>HEC</states>\n" +
                                                 "\t\t<observations>ZWD" +
                                                 "</observations" +
                                                 ">\n" +
                                                 "\t\t<order>5</order>\n" +
                                                 "\t</meta>\n" +
                                                 "</markovmodel>\n";

    private static final String xmlDefinitionTrained = "<markovmodel>\n" +
                                                        "\t<meta>\n" +
                                                        "\t\t<states>HEC</states>\n" +
                                                        "\t\t<observations>ZWD" +
                                                        "</observations" +
                                                        ">\n" +
                                                        "\t\t<order>5</order>\n" +
                                                        "\t</meta>\n" +
                                                        "\t<ortho>\n" +
                                                        "\t</ortho>\n" +
                                                        "</markovmodel>\n";


    /**
     * Parse the parameters of an markovmodel.xml file for an HMM object as HMMParams object.
     *
     * @param hmmFile the string to the markovmodel file
     * @return a HMMParams object containing all relevant parameters for the HMM
     */
    public static HMMParams parseXML(java.lang.String hmmFile)
    {
        HMMParams params = HMMParams.newInstance();
        SAXBuilder builder = new SAXBuilder();
        Document document;
        try
        {
            _LOGGER.info("Parsing markovmodel xml.");
            document = builder.build(hmmFile);
            setStandardParams(document, params);
            setTrainingParams(document, params);
        }
        catch (IOException | JDOMException e)
        {
            _LOGGER.error("Could not open file: " + e.getMessage());
        }
        return params;
    }

    @SuppressWarnings("unchecked")
    private static void setTrainingParams(Document document, HMMParams params)
    {
        _LOGGER.info("Parsing training parameters.");
        String exit = "Your XML format is wrong! It should look like this:\n" + xmlDefinitionTrained;
        Element rootNode = document.getRootElement();
        Element ortho = rootNode.getChild("ortho");
        if (ortho == null)
            return;
        if (ortho.getChild("starts") == null ||
            ortho.getChild("emissions") == null ||
            ortho.getChild("transitions") == null)
        {
            _LOGGER.error("Some elements in the xml are null.");
            net.digital_alexandria.lvm4j.util.System.exit(exit, ExitCode.EXIT_ERROR);
        }

        Element start = ortho.getChild("starts");
        List list = start.getChildren();
        for (int i = 0; i < list.size(); i++)
        {
            Element node = (Element) list.get(i);
            String state = node.getAttribute("state").getValue();
            double prob = Double.parseDouble(node.getText());
            params.startProbabilities().add(new Pair<>(state, prob));
        }

        Element transitions = ortho.getChild("transitions");
        list = transitions.getChildren();
        for (int i = 0; i < list.size(); i++)
        {
            Element node = (Element) list.get(i);
            String source = node.getAttribute("source").getValue();
            String sink = node.getAttribute("sink").getValue();
            double prob = Double.parseDouble(node.getText());
            params.transitionProbabilities().add(new Triple<>(source, sink, prob));
        }

        Element emissions = ortho.getChild("emissions");
        list = emissions.getChildren();
        for (int i = 0; i < list.size(); i++)
        {
            Element node = (Element) list.get(i);
            String source = node.getAttribute("source").getValue();
            String sink = node.getAttribute("sink").getValue();
            double prob = Double.parseDouble(node.getText());
            params.emissionProbabilities().add(new Triple<>(source, sink, prob));
        }
        params.setTrainingParam(true);
    }

    private static void setStandardParams(Document document, HMMParams params)
    {
        _LOGGER.info("Parsing standard parameters.");
        String exit = "Your XML format is wrong! It should look like this:\n" + xmlDefinition;
        Element rootNode = document.getRootElement();
        Element meta = rootNode.getChild("meta");
        if (meta == null ||
            meta.getChild("states") == null ||
            meta.getChild("observations") == null ||
            meta.getChild("order") == null)
                net.digital_alexandria.lvm4j.util.System.exit(exit,  ExitCode.EXIT_ERROR);
        char states[] = meta.getChild("states").getValue().toCharArray();
        char observations[] = meta.getChild("observations").getValue()
                                  .toCharArray();
        int order = Integer.parseInt(meta.getChild("order").getValue());
        if (states.length == 0 || observations.length == 0)
            net.digital_alexandria.lvm4j.util.System.exit(exit,  ExitCode.EXIT_ERROR);
        params.observations(observations);
        params.order(order);
        params.states(states);
    }

    /**
     * Write the HMM parameters to a xml file.
     *
     * @param hmm the markovmodel of which should be written
     * @param file  the output file
     */
    public static void writeXML(HMM hmm, String file)
    {
        try
        {
            _LOGGER.info("Writing markovmodel to xml.");
            Element hmmxml = new Element("markovmodel");
            Document doc = new Document(hmmxml);
            doc.setRootElement(hmmxml);

            Element meta = new Element("meta");
            hmmxml.addContent(meta);
            addMeta(hmm, meta);
            if (hmm.isTrained())
            {
                Element ortho = new Element("ortho");
                hmmxml.addContent(ortho);
                addOrtho(hmm, ortho);
            }
            XMLOutputter xmlOutput = new XMLOutputter();
            xmlOutput.setFormat(Format.getPrettyFormat());
            xmlOutput.output(doc, new FileWriter(file));
        }
        catch (IOException e)
        {
            _LOGGER.error("Error when writing xmlFile: " + e.getMessage());
        }
    }

    @SuppressWarnings("unchecked")
    private static void addOrtho(HMM ssHMM, Element ortho)
    {
        _LOGGER.info("Writing ortho information (trained parameters).");
        Element starts = new Element("starts");
        ortho.addContent(starts);
        ssHMM.states().stream().filter(s -> s.state().length() == 1).forEach(s -> {
            Element start = new Element("start");
            starts.addContent(start);
            start.setAttribute("state", s.state());
            start.setText(String.valueOf(s.startingProbability()));
        });
        Element transitions = new Element("transitions");
        ortho.addContent(transitions);
        for (WeightedArc t : ssHMM.transitions())
            add(t, "transition", transitions);
        Element emissions = new Element("emissions");
        ortho.addContent(emissions);
        for (WeightedArc e : ssHMM.emissions())
            add(e, "emission", emissions);
    }

    @SuppressWarnings("unchecked")
    private static void add(WeightedArc arc, String lab, Element elem)
    {
        HMMNode<Character, String> src = (HMMNode<Character, String>) arc.source();
        HMMNode<Character, String> sink = (HMMNode<Character, String>) arc.sink();
        Element el = new Element(lab);
        elem.addContent(el);
        el.setAttribute("source", src.state());
        el.setAttribute("sink", sink.state());
        el.setText(String.valueOf(arc.weight()));
    }

    private static void addMeta(HMM ssHMM, Element meta)
    {
        _LOGGER.info("Writing meta information (trained parameters).");
        Set<String> sb = ssHMM.states().stream()
                               .map(s -> String.valueOf(s.label()))
                               .collect(Collectors.toSet());
        StringBuilder states = new StringBuilder();
        sb.stream().forEach(states::append);

        sb = ssHMM.observations().stream()
                  .map(s -> String.valueOf(s.label()))
                  .collect(Collectors.toSet());
        StringBuilder observations = new StringBuilder();
        sb.forEach(observations::append);
        meta.addContent(
            new Element("states").setText(states.toString()));
        meta.addContent(
            new Element("observations").setText(observations.toString()));
        meta.addContent(
            new Element("order").setText(String.valueOf(ssHMM.order())));
    }
}
