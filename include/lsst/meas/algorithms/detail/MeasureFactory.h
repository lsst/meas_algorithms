// -*- LSST-C++ -*-
#if !defined(LSST_DETECTION_MEASURE_FACTORY_H)
#define LSST_DETECTION_MEASURE_FACTORY_H

#include "boost/cstdint.hpp"
#include "boost/shared_ptr.hpp"
#include "lsst/pex/logging/Log.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/ImageUtils.h"

namespace lsst {
namespace meas {
namespace algorithms {

/************************************************************************************************************/
/**
 * Return a MeasureProperty of the requested variety
 *
 * Must be here, as it's declared a friend by MeasureProperty
 *
 * @throws std::runtime_error if name can't be found
 */
template<typename MeasurePropertyT, typename ImageT>
MeasurePropertyT* createMeasureProperty(
        std::string const& name,        ///< desired variety
        boost::shared_ptr<ImageT const> image=boost::shared_ptr<ImageT const>(), ///< the image to process
        MeasurePropertyT const* =NULL   ///< a MeasurePropertyT to disambiguate the function
                                              )
{
    return MeasurePropertyT::lookup(name).create(image);
}

/************************************************************************************************************/
/**
 * A polymorphic base class for MeasureProperty factories
 */
template<typename MeasurePropertyT>
class MeasurePropertyFactoryBase : public lsst::daf::base::Citizen {
public:
    typedef MeasurePropertyT MeasureProperty;
    MeasurePropertyFactoryBase() : lsst::daf::base::Citizen(typeid(this)) {}
    virtual ~MeasurePropertyFactoryBase() {}
    virtual MeasurePropertyT *
    create(typename MeasurePropertyT::ImageT::ConstPtr=typename MeasurePropertyT::ImageT::ConstPtr()) = 0;
};

/**
 * Create a particular sort of MeasureProperty
 */
#if !defined(SWIG)
template<template<typename T> class MeasurePropertyT, typename ImageT>
class MeasurePropertyFactory
    : public MeasurePropertyFactoryBase<typename MeasurePropertyT<ImageT>::MeasurePropertyBase >
{
public:
    typedef MeasurePropertyFactoryBase<MeasurePropertyT<ImageT> > baseClass;
    /**
     * Return a new MeasurePropertyT<ImageT>
     */
    MeasurePropertyT<ImageT> *create(typename ImageT::ConstPtr image) {
        return new MeasurePropertyT<ImageT>(image);
    }
};
#endif

/************************************************************************************************************/
/**
 * @brief A pure virtual base class to calculate some property of an image
 *
 * Different implementations will use different algorithms
 */
template<typename T, typename _ImageT>
class MeasureProperty : public boost::noncopyable {
public:
    typedef _ImageT ImageT;

    MeasureProperty(boost::shared_ptr<ImageT const> image=boost::shared_ptr<ImageT const>()) : _image(image)
    {
        static bool _registered = false;

        if (!_registered) {
#if 0                                   // We don't actually declare the (pure virtual) base class
            MeasureProperty::declare("base", new MeasurePropertyFactory<MeasureProperty>());
#endif
            _registered = true;
        }        
    }
    virtual ~MeasureProperty() {}

    /**
     * Return the image that we're working on
     */
    boost::shared_ptr<ImageT const> getImage() const {
        return _image;
    }
    /**
     * Change the image that we're working on
     */
    void setImage(boost::shared_ptr<ImageT const> image) const {
        _image = image;
    }
protected:
#if !defined(SWIG)
    friend T *createMeasureProperty<T, ImageT>(std::string const& name, boost::shared_ptr<ImageT const> image, T const*);
#endif

    /**
     * Declare a MeasurePropertyFactory for a variety "name"
     *
     * @throws std::runtime_error if name is already declared
     */
    static void declare(std::string name,                        ///< name of variety
                        MeasurePropertyFactoryBase<T>* factory ///< Factory to make this sort of widget
                       ) {
        (void)_registry(name, factory);
    }

    /**
     * Return the named MeasurePropertyFactory
     *
     * @throws std::runtime_error if name can't be found
     */
    static MeasurePropertyFactoryBase<T>& lookup(std::string name ///< desired variety
                                                     ) {
        return _registry(name, NULL);
    }

    /**
     * Register the factory that builds a particular sort of MeasureProperty
     *
     * \note This function returns bool so that it can be used in an initialisation at file scope to do the
     * actual registration
     */
#if !defined(SWIG)
    template<template<typename TT> class MeasurePropertyT, typename ImageT>
    friend bool registerMe(std::string const& name);
#endif
private:
    mutable boost::shared_ptr<ImageT const> _image;
    static MeasurePropertyFactoryBase<T>& _registry(std::string name,
                                                    MeasurePropertyFactoryBase<T>* factory);
};

/************************************************************************************************************/

#if !defined(SWIG)
template<template<typename T> class MeasurePropertyT, typename ImageT>
bool registerMe(std::string const& name) {
    static bool _registered = false;
    
    //std::cout << "Registering " << name.c_str() << " " << typeid(MeasurePropertyT<ImageT>).name() << " " <<  _registered << std::endl;

    if (!_registered) {
        typedef MeasurePropertyFactory<MeasurePropertyT, ImageT> Factory;

        Factory *factory = new Factory();
        factory->markPersistent();

        MeasurePropertyT<ImageT>::declare(name, factory);
        _registered = true;
    }
    
    return true;
}
#endif

/************************************************************************************************************/
/*
 * Register a factory object by name;  if the factory's NULL, return the named factory
 */
template<typename T, typename ImageT>
MeasurePropertyFactoryBase<T>&
MeasureProperty<T, ImageT>::_registry(std::string name,
                                      MeasurePropertyFactoryBase<T>* factory)
{
    static std::map<std::string const, MeasurePropertyFactoryBase<T> *> _theRegistry;

    //std::cout << "Looking up  " << name.c_str() << " " << typeid(T).name() << " theReg: " << &_theRegistry << " " <<  _theRegistry.size() << std::endl;

    typename std::map<std::string const, MeasurePropertyFactoryBase<T> *>::iterator el = _theRegistry.find(name);    
    if (el == _theRegistry.end()) {        // failed to find name
        if (factory) {
            _theRegistry[name] = factory;
        } else {
            throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundException, 
                              "MeasureProperty of type \"" + name + "\" is not implemented");
        }
    } else {
        if (!factory) {
            factory = (*el).second;
        } else if(factory == (*el).second) {
            ;                           // OK
        } else {
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, 
                              "MeasureProperty of type \"" + name + "\" is already declared");
        }
    }

    return *factory;
}

}}}
#endif
