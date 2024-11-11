#ifndef IMAGE_IO_HPP
#define IMAGE_IO_HPP
#include <algorithm>
#include <boost/gil.hpp>
#include <boost/gil/extension/dynamic_image/any_image.hpp>
#include <boost/gil/extension/dynamic_image/apply_operation.hpp>
#include <boost/gil/extension/io/bmp/tags.hpp>
#include <boost/gil/extension/io/jpeg/tags.hpp>
#include <boost/gil/io/io.hpp>
#include <boost/gil/extension/io/jpeg.hpp>
#include <boost/gil/extension/io/png.hpp>
#include <boost/gil/extension/io/bmp.hpp>
#include <boost/gil/io/write_view.hpp>
#include <boost/gil/typedefs.hpp>
#include <boost/variant2/variant.hpp>

template <typename Image>
auto gray_image_from_file(std::string filename)
{
    using AnyImage = boost::gil::any_image<
    boost::gil::gray8_image_t,
    boost::gil::rgb8_image_t,
    boost::gil::rgba8_image_t
    >;

    AnyImage image;

    // Определите формат по расширению или по сигнатуре
    if ((filename.find(".jpg") != filename.npos) || (filename.find(".jpeg") != filename.npos)) {
        read_image(filename, image, boost::gil::jpeg_tag());
    } else if (filename.find(".png") != filename.npos) {
        read_image(filename, image, boost::gil::png_tag());
    } else if (filename.find(".bmp") != filename.npos) {
        read_image(filename, image, boost::gil::bmp_tag());
    } else {
        throw std::runtime_error("Unsupported image format.");
    }

    Image out(image.width(),image.height());
    
    boost::variant2::visit
    ([&](const auto& img_view) {
        using ViewType = std::decay_t<decltype(img_view)>;
        
        if constexpr (std::is_same_v<ViewType, boost::gil::gray8_view_t>) {
            std::cout << "Image format: Grayscale 8-bit\n";
            for (auto i =0; i < img_view.width(); ++i) 
            for (auto j = 0; j < img_view.height(); ++j)
            out(i,j) = img_view(i,j);
        } else if constexpr (std::is_same_v<ViewType, boost::gil::rgb8_view_t>) {
            std::cout << "Image format: RGB 8-bit\n";
            for (auto i =0; i < img_view.width(); ++i) 
            for (auto j = 0; j < img_view.height(); ++j)
            out(i,j) = 0.299*img_view(i,j)[0] + 0.587*img_view(i,j)[1] + 0.114*img_view(i,j)[2];
        } else if constexpr (std::is_same_v<ViewType, boost::gil::rgba8_view_t>) {
            std::cout << "Image format: RGBA 8-bit\n";
            for (auto i =0; i < img_view.width(); ++i) 
            for (auto j = 0; j < img_view.height(); ++j)
            out(i,j) = 0.299*img_view(i,j)[0] + 0.587*img_view(i,j)[1] + 0.114*img_view(i,j)[2];
        } else {
            std::cout << "Unknown format\n";
        }
    },(view(image)));

    return out;
}


template <typename Image>
auto gray_image_to_file(Image im, unsigned w, unsigned h, std::string filename)
{
    boost::gil::gray8_image_t img(w,h);

    for (auto i =0u; i < w; ++i)
    {
        for (auto j = 0u; j < h; ++j)
        {
            img._view(i,j) = im(i,j);
        }
    }

    // Определите формат по расширению или по сигнатуре
    if ((filename.find(".jpg") != filename.npos) || (filename.find(".jpeg") != filename.npos)) {
        boost::gil::write_view(filename, img._view, boost::gil::jpeg_tag());
    } else if (filename.find(".png") != filename.npos) {
        boost::gil::write_view(filename, img._view, boost::gil::png_tag());
    } else if (filename.find(".bmp") != filename.npos) {
        // boost::gil::write_view(filename, img._view, boost::gil::bmp_tag());
        throw std::runtime_error("Unsupported image format.");
    } else {
        throw std::runtime_error("Unsupported image format.");
    }
}


#endif // IMAGE_IO_HPP
