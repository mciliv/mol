import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import 'package:image_picker/image_picker.dart';
import 'package:flutter_svg/flutter_svg.dart';
import 'dart:io';
import 'payment_service.dart';
import 'analysis_service.dart';

class PhotoPicker extends StatefulWidget {
  @override
  _PhotoPickerState createState() => _PhotoPickerState();
}

class _PhotoPickerState extends State<PhotoPicker> {
  final _urlController = TextEditingController();
  final _imagePicker = ImagePicker();
  
  @override
  Widget build(BuildContext context) {
    return Consumer2<PaymentService, AnalysisService>(
      builder: (context, paymentService, analysisService, child) {
        return Padding(
          padding: EdgeInsets.symmetric(vertical: 20),
          // FLOWING STRUCTURE - Column instead of rigid Row
          child: Column(
            mainAxisSize: MainAxisSize.min, // KEY: Prevents overflow
            children: [
              // Upload button - clean, no box
              _buildUploadButton(paymentService, analysisService),
              
              SizedBox(height: 30),
              
              // URL input section - flowing
              _buildUrlSection(paymentService, analysisService),
            ],
          ),
        );
      },
    );
  }
  
  Widget _buildUploadButton(PaymentService paymentService, AnalysisService analysisService) {
    return GestureDetector(
      onTap: () => _pickImage(paymentService, analysisService),
      child: Padding(
        padding: EdgeInsets.all(20),
        child: Column(
          mainAxisSize: MainAxisSize.min, // KEY: Prevents overflow
          children: [
            SvgPicture.asset(
              'assets/photo-upload.svg',
              width: 24,
              height: 24,
              colorFilter: ColorFilter.mode(
                Colors.white.withOpacity(0.7),
                BlendMode.srcIn,
              ),
            ),
            SizedBox(height: 12),
            Text(
              'Upload',
              style: TextStyle(
                color: Colors.white.withOpacity(0.7),
                fontSize: 16,
                fontWeight: FontWeight.w500,
              ),
            ),
          ],
        ),
      ),
    );
  }
  
  Widget _buildUrlSection(PaymentService paymentService, AnalysisService analysisService) {
    return Column(
      mainAxisSize: MainAxisSize.min, // KEY: Prevents overflow
      children: [
        Text(
          'Or paste image URL:',
          style: TextStyle(
            color: Colors.white.withOpacity(0.7),
            fontSize: 14,
          ),
        ),
        SizedBox(height: 12),
        Row(
          children: [
            Expanded(
              child: TextField(
                controller: _urlController,
                decoration: InputDecoration(
                  hintText: 'Image URL...',
                  isDense: true,
                  contentPadding: EdgeInsets.symmetric(horizontal: 16, vertical: 12),
                ),
                style: TextStyle(color: Colors.white, fontSize: 14),
                onSubmitted: (_) => _analyzeUrl(paymentService, analysisService),
              ),
            ),
            SizedBox(width: 12),
            GestureDetector(
              onTap: () => _analyzeUrl(paymentService, analysisService),
              child: Padding(
                padding: EdgeInsets.all(12),
                child: Icon(
                  Icons.search,
                  size: 16,
                  color: Colors.white.withOpacity(0.7),
                ),
              ),
            ),
          ],
        ),
      ],
    );
  }
  
  Future<void> _pickImage(PaymentService paymentService, AnalysisService analysisService) async {
    if (!paymentService.canAnalyze()) {
      paymentService.showPayment();
      return;
    }
    
    try {
      final XFile? image = await _imagePicker.pickImage(
        source: ImageSource.gallery,
        maxWidth: 1920,
        maxHeight: 1080,
        imageQuality: 85,
      );
      
      if (image != null) {
        final imageFile = File(image.path);
        await analysisService.analyzeImage(imageFile);
        await paymentService.incrementUsage();
        
        ScaffoldMessenger.of(context).showSnackBar(
          SnackBar(content: Text('Image analysis complete')),
        );
      }
    } catch (error) {
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('Image pickup failed: $error')),
      );
    }
  }
  
  Future<void> _analyzeUrl(PaymentService paymentService, AnalysisService analysisService) async {
    final url = _urlController.text.trim();
    if (url.isEmpty) return;
    
    if (!paymentService.canAnalyze()) {
      paymentService.showPayment();
      return;
    }
    
    try {
      _urlController.clear();
      await analysisService.analyzeImageUrl(url);
      await paymentService.incrementUsage();
      
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('URL analysis complete')),
      );
    } catch (error) {
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('URL analysis failed: $error')),
      );
    }
  }
  
  @override
  void dispose() {
    _urlController.dispose();
    super.dispose();
  }
} 