package edu.cmu.tetradapp.editor.factory;


import edu.cmu.tetrad.annotation.AlgorithmDescription;
import edu.cmu.tetradapp.editor.AlgorithmDescriptionClass;
import org.reflections.Reflections;

import java.lang.annotation.Annotation;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * Author : Jeremy Espino MD
 * Created  6/30/17 11:20 AM
 */
public class AlgorithmDescriptionFactory {
    private static AlgorithmDescriptionFactory ourInstance = new AlgorithmDescriptionFactory();

    public static AlgorithmDescriptionFactory getInstance() {
        return ourInstance;
    }

    private AlgorithmDescriptionFactory() {
    }

    public List<AlgorithmDescriptionClass> getAlgorithmDescriptions() {

        ArrayList<AlgorithmDescriptionClass> algorithmDescriptions = new ArrayList<AlgorithmDescriptionClass>();

        // find all classes that implement algorithm
        Reflections reflections = new Reflections("edu.cmu.tetrad.algcomparison");
        Set<Class<?>> classes = reflections.getTypesAnnotatedWith(edu.cmu.tetrad.annotation.AlgorithmDescription.class);

        for (Class clazz : classes) {
            AlgorithmDescriptionClass algorithmDescription;
            Annotation annotation = clazz.getAnnotation(AlgorithmDescription.class);

            if (annotation instanceof AlgorithmDescription) {
                AlgorithmDescription myAnnotation = (AlgorithmDescription) annotation;

                algorithmDescription = new AlgorithmDescriptionClass(myAnnotation.name(), myAnnotation.algType(), myAnnotation.oracleType());
                algorithmDescriptions.add(algorithmDescription);
            }
        }


        return algorithmDescriptions;
    }

    public static void main(String[] args) {
        AlgorithmDescriptionFactory.getInstance().getAlgorithmDescriptions();
    }
}
